async function load_grid_from_json() {
  const response = await fetch("bar_density_grid.json");
  if (!response.ok) throw new Error("Failed to load bar density grid");
  const data = await response.json();
  console.log("Bar density grid loaded.");
  return data;
}

/**
 * Calculates the overlap ratio between two boundary boxes.
 * @param {object} cellBounds - The bounds of the grid cell { south, west, north, east }.
 * @param {object} gameBounds - The bounds of the custom game area { south, west, north, east }.
 * @returns {number} The overlap ratio, from 0.0 (no overlap) to 1.0 (full overlap).
 */
function calculateOverlapRatio(cellBounds, gameBounds) {
    // 1. Find the boundaries of the intersection rectangle
    const intersectSouth = Math.max(cellBounds.south, gameBounds.south);
    const intersectWest = Math.max(cellBounds.west, gameBounds.west);
    const intersectNorth = Math.min(cellBounds.north, gameBounds.north);
    const intersectEast = Math.min(cellBounds.east, gameBounds.east);

    // 2. Check if there is any overlap
    // If north < south or east < west, the rectangles do not overlap.
    if (intersectNorth <= intersectSouth || intersectEast <= intersectWest) {
        return 0.0;
    }

    // 3. Calculate areas in "degree-squared" units
    // The distortion (cos(lat)) cancels out when we take the ratio.
    const intersectArea = (intersectNorth - intersectSouth) * (intersectEast - intersectWest);
    const cellArea = (cellBounds.north - cellBounds.south) * (cellBounds.east - cellBounds.west);

    // 4. Return the ratio
    // (Handle potential division by zero, though cellArea should always be > 0)
    return cellArea > 0 ? intersectArea / cellArea : 0.0;
}

/**
 * Generates a list of weighted cells based on their bar count
 * and partial overlap with the game bounds.
 * * @param {object} barDensityGrid - The full { 'lat_lon': count } object.
 * @param {object} gameBounds - The custom game area { south, west, north, east }.
 * @param {number} gridStep - The size of your grid cells (e.g., 1.0 for 1x1).
 * @returns {object} An object containing { weightedCells, totalWeight }.
 */
function getAdjustedWeightedCells(barDensityGrid, gameBounds, gridStep) {
    const weightedCells = [];
    let totalWeight = 0.0;

    for (const key in barDensityGrid) {
        const count = barDensityGrid[key];

        // Skip empty cells immediately
        if (count === 0) {
            continue;
        }

        // Re-create the cell's bounds from its key
        const [lat_s, lon_w] = key.split('_').map(parseFloat);
        const cellBounds = {
            south: lat_s,
            west: lon_w,
            north: lat_s + gridStep,
            east: lon_w + gridStep
        };

        // Calculate the overlap
        const overlapRatio = calculateOverlapRatio(cellBounds, gameBounds);

        // If the cell overlaps, calculate its new weight and add it
        if (overlapRatio > 0.0) {
            // This is the crucial part:
            const adjustedWeight = count * overlapRatio;

            weightedCells.push({
                key: key,
                bounds: cellBounds,
                originalCount: count,
                weight: adjustedWeight // Use this new weight for selection
            });
            
            totalWeight += adjustedWeight;
        }
    }
    
    // We round the total weight slightly to avoid floating point issues
    return { weightedCells, totalWeight: parseFloat(totalWeight.toFixed(5)) };
}

/**
 * Picks a random cell from the list, respecting the new adjusted weights.
 * @param {Array} weightedCells - The list from getAdjustedWeightedCells.
 *V @param {number} totalWeight - The total weight from getAdjustedWeightedCells.
 * @returns {object} The chosen cell object.
 */
function pickRandomCell(weightedCells, totalWeight) {
    if (weightedCells.length === 0 || totalWeight <= 0) {
        throw new Error("No weighted cells found in this area.");
    }

    let randomNum = Math.random() * totalWeight;

    for (const cell of weightedCells) {
        // Use the new decimal weight
        if (randomNum < cell.weight) {
            return cell; // This is the chosen cell
        }
        randomNum -= cell.weight;
    }

    // Fallback for floating point errors
    return weightedCells[weightedCells.length - 1];
}

// --- Your new logic ---

export async function drawRandomBarFromDensityGrid(GAME_BOUNDS) {
  const GRID_STEP = 2.0;

  const barDensityGrid = await load_grid_from_json(); // could be loaded in main file!!

  const { weightedCells, totalWeight } = getAdjustedWeightedCells(
    barDensityGrid,
    GAME_BOUNDS,
    GRID_STEP
  );

  // 1. Get weighted cells inside your large GAME_BOUNDS
  if (weightedCells.length === 0) {
    throw new Error("No bars found in any grid cell in this area.");
  }

  // 2. Pick one cell, weighted by its bar count
  const chosenCell = pickRandomCell(weightedCells, totalWeight);

  // 3. IMPORTANT: Use the chosen cell's bounds for the Overpass query.
  // We must also clip these bounds to be *within* the original GAME_BOUNDS
  // in case the cell was only partially overlapping.
  const queryBounds = {
    south: Math.max(chosenCell.bounds.south, GAME_BOUNDS.south),
    west: Math.max(chosenCell.bounds.west, GAME_BOUNDS.west),
    north: Math.min(chosenCell.bounds.north, GAME_BOUNDS.north),
    east: Math.min(chosenCell.bounds.east, GAME_BOUNDS.east),
  };

  // 4. Call your *original* fetch function, but with the new, smaller bounds
  // I've modified it to accept bounds as an argument:
  const barsInCell = await fetchBarsInBounds(queryBounds);

  return barsInCell;
}

// --- End new logic ---

async function fetchBarsInBounds(bounds) {
  const overpassUrl = "https://overpass-api.de/api/interpreter";
  const query = `
        [out:json][timeout:100];
        (
            node["amenity"="bar"](${bounds.south},${bounds.west},${bounds.north},${bounds.east});
            node["amenity"="pub"](${bounds.south},${bounds.west},${bounds.north},${bounds.east});
        );
        out body;
    `;

  const response = await fetch(overpassUrl, {
    method: "POST",
    body: query,
  });

  if (!response.ok) throw new Error("Failed to fetch bars from OpenStreetMap");
  const data = await response.json();

  const bars = data.elements
    .filter((el) => el.tags && el.tags.name)
    .map((el) => ({
      name: el.tags.name,
      lat: el.lat,
      lng: el.lon,
      inside: true,
    }));

  if (bars.length === 0) {
    // If the chosen cell has 0 named bars (even if the grid said it had 10),
    // you might need to retry by picking another cell.
    // For now, we'll just throw an error.
    throw new Error("No named bars found in the randomly selected area.");
  }

  return bars;
}
