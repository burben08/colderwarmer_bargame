async function load_grid_from_json() {
  const response = await fetch("bar_density_grid.json");
  if (!response.ok) throw new Error("Failed to load bar density grid");
  const data = await response.json();
  console.log("Bar density grid loaded.");
  return data;
}

async function getWeightedCells(GAME_BOUNDS) {
  const grid = await load_grid_from_json();
  const weightedCells = [];
  let totalWeight = 0;

  console.log(grid);

  // Loop through the grid
  for (const key in grid) {
    const [lat_s, lon_w] = key.split("_").map(parseFloat);
    const lat_n = lat_s + 2.0; // Assuming 2.0 step
    const lon_e = lon_w + 2.0; // Assuming 2.0 step

    // Check if this cell *overlaps* with the game bounds
    const overlaps =
      lat_s < GAME_BOUNDS.north &&
      lat_n > GAME_BOUNDS.south &&
      lon_w < GAME_BOUNDS.east &&
      lon_e > GAME_BOUNDS.west;

    const count = grid[key];

    if (overlaps && count > 0) {
      weightedCells.push({
        key: key,
        bounds: { south: lat_s, west: lon_w, north: lat_n, east: lon_e },
        count: count,
      });
      totalWeight += count;
    }
  }
  return { weightedCells, totalWeight };
}

function pickRandomCell(weightedCells, totalWeight) {
  let randomNum = Math.random() * totalWeight;

  for (const cell of weightedCells) {
    if (randomNum < cell.count) {
      return cell; // This is the chosen cell
    }
    randomNum -= cell.count;
  }
  // Fallback in case of rounding errors
  return weightedCells[weightedCells.length - 1];
}

// --- Your new logic ---

export async function drawRandomBarFromDensityGrid(GAME_BOUNDS) {
  // 1. Get weighted cells inside your large GAME_BOUNDS
  const { weightedCells, totalWeight } = await getWeightedCells(GAME_BOUNDS);
  console.log("Weighted cells found:", weightedCells.length);
  console.log("Total weight:", totalWeight);
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

  // 5. Pick a random bar from this much smaller list
  const randomBar = barsInCell[Math.floor(Math.random() * barsInCell.length)];

  return barsInCell;
}

// --- End new logic ---

async function fetchBarsInBounds(bounds) {
  const overpassUrl = "https://overpass-api.de/api/interpreter";
  const query = `
        [out:json][timeout:25];
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
