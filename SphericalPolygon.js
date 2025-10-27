// --- Module-level private variables ---
// These will be set by the exported buildPolygon function
let CURRENT_GAME_BOUNDS = {};
let CURRENT_GAME_BOUNDS_CORNERS_LONLAT = [];

// --- Helper math ---
const toRad = d => d * Math.PI / 180
const toDeg = r => r * 180 / Math.PI

// Input: lon, lat in degrees
// Output: [x, y, z] cartesian vector
export function lonLatToCartesian(lon, lat) {
  const φ = toRad(lat),
    λ = toRad(lon)
  return [Math.cos(φ) * Math.cos(λ), Math.cos(φ) * Math.sin(λ), Math.sin(φ)]
}
// Input: [x, y, z] cartesian vector
// Output: [lon, lat] in degrees
function cartesianToLonLat(v) {
  const [x, y, z] = v
  const lon = Math.atan2(y, x)
  const hyp = Math.sqrt(x * x + y * y)
  const lat = Math.atan2(z, hyp)
  return [toDeg(lon), toDeg(lat)]
}

// --- Vector math ---
const cross = (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
const norm = v => {
  const m = Math.hypot(v[0], v[1], v[2]);
  // Avoid division by zero
  if (m < 1e-15) return [0, 0, 0];
  return v.map(x => x / m)
}
function dot(a, b) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
function sub(a, b) {
  return [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
function add(a, b) {
  return [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}
function mul(a, s) {
  return [a[0] * s, a[1] * s, a[2] * s]
}
function length(a) {
  return Math.hypot(a[0], a[1], a[2])
}
function normalizeLon(l) {
  let x = ((l + 180) % 360 + 360) % 360 - 180;
  return x
}
// Input: p1, p2 as [x, y, z] cartesian vectors
export function sphericalDistance(p1, p2) {
  const r = length(p1) // Assumes points are on a sphere of radius r
  const r2 = length(p2) // Should be same radius
  const cosTheta = dot(p1, p2) / (r * r2)
  const theta = Math.acos(Math.max(-1, Math.min(1, cosTheta)))
  return r * theta // Returns arc distance
}

// Input: Two lines, L1 and L2, each defined by two points in [lon, lat] format
// Output: Intersection point in [lon, lat] in degrees
function getIntersectionOnSphere(L1, L2) {
  const A1 = lonLatToCartesian(L1[0][0], L1[0][1]);
  const B1 = lonLatToCartesian(L1[1][0], L1[1][1]);
  const A2 = lonLatToCartesian(L2[0][0], L2[0][1]);
  const B2 = lonLatToCartesian(L2[1][0], L2[1][1]);
  const n1 = norm(cross(A1, B1))
  const n2 = norm(cross(A2, B2))
  const p = norm(cross(n1, n2)) // Intersection point (or its antipode)
  const φ = Math.asin(p[2])
  const λ = Math.atan2(p[1], p[0])
  return [toDeg(λ), toDeg(φ)]
}

// Input: p1LonLat, p2LonLat as [lon, lat]
// Output: Two points [M, Q] defining the bisection great circle, each as [lon, lat]
function getBisectionLine(p1LonLat, p2LonLat) {
  const A = lonLatToCartesian(p1LonLat[0], p1LonLat[1]);
  const B = lonLatToCartesian(p2LonLat[0], p2LonLat[1]);
  const n = cross(A, B); // Normal to the p1-p2 great circle

  const m = norm([A[0] + B[0], A[1] + B[1], A[2] + B[2]]); // Midpoint vector
  const t = cross(n, m); // Tangent at midpoint, on the p1-p2 great circle
  const q = norm(n); // Normal to p1-p2 great circle, a point on the bisector

  const M = cartesianToLonLat(m); // Midpoint in [lon, lat]
  const Q = cartesianToLonLat(q); // Point 'q' in [lon, lat]

  return [M, Q]; // M and Q define the bisector great circle
}

export function getLineCoordsBetween(p1LonLat, p2LonLat, steps = 300) {
  let [lon1, lat1] = p1LonLat
  let [lon2, lat2] = p2LonLat

  const A = lonLatToCartesian(lon1, lat1)
  const B = lonLatToCartesian(lon2, lat2)
  const coords = []

  for (let i = 0; i <= steps; i++) {
    const t = i / steps
    // Linear interpolation in 3D cartesian space
    const v = norm([
      A[0] * (1 - t) + B[0] * t,
      A[1] * (1 - t) + B[1] * t,
      A[2] * (1 - t) + B[2] * t
    ])
    let [lon, lat] = cartesianToLonLat(v)
    coords.push([lon, lat])
  }

  return coords
}

// Get intersections with true sphere boundaries
// Input: p1LonLat, p2LonLat as [lon, lat]
// Output: Array of intersection points as [lon, lat]
function getBoundaryIntersections(p1LonLat, p2LonLat) {
  const geodesicBisection1 = getBisectionLine(p1LonLat, p2LonLat);
  const geodesicBisection2 = getBisectionLine(p2LonLat, p1LonLat); // Antipodal points for bisection line

  // Boundary lines defined by [lon, lat] points
  // Reads from the module-level variable
  console.log('CURRENT_GAME_BOUNDS_CORNERS_LONLAT:', CURRENT_GAME_BOUNDS_CORNERS_LONLAT);
  const boundaryLines = [
    [CURRENT_GAME_BOUNDS_CORNERS_LONLAT[1], CURRENT_GAME_BOUNDS_CORNERS_LONLAT[0]], // West
    [CURRENT_GAME_BOUNDS_CORNERS_LONLAT[2], CURRENT_GAME_BOUNDS_CORNERS_LONLAT[1]], // North
    [CURRENT_GAME_BOUNDS_CORNERS_LONLAT[3], CURRENT_GAME_BOUNDS_CORNERS_LONLAT[2]], // East
    [CURRENT_GAME_BOUNDS_CORNERS_LONLAT[0], CURRENT_GAME_BOUNDS_CORNERS_LONLAT[3]]  // South
  ];

  const intersections = [];

  for (const line of boundaryLines) {
    const i1 = getIntersectionOnSphere(geodesicBisection1, line)
    const i2 = getIntersectionOnSphere(geodesicBisection2, line)
    intersections.push(i1, i2)
  }

  return intersections;
}

// Input: p1LonLat, p2LonLat as [lon, lat], latC = constant latitude in degrees
// Output: Array of longitudes [lon1, lon2] where geodesic intersects latC, or null
function getIntersectionWithConstantLat(p1LonLat, p2LonLat, latC) {
  // Convert to Cartesian
  const pA = lonLatToCartesian(p1LonLat[0], p1LonLat[1])
  const pB = lonLatToCartesian(p2LonLat[0], p2LonLat[1])

  // Get great circle normal
  const n = cross(pA, pB)
  const nLen = length(n)

  if (nLen < 1e-12) {
    // Points are antipodal or identical
    return null
  }
  const nNorm = mul(n, 1 / nLen)

  // We want points where n·p = 0 AND z = sin(latC)
  // Build orthonormal basis (u, v) in the great circle plane
  let u = pA
  u = sub(u, mul(nNorm, dot(u, nNorm))) // Make u perpendicular to n
  let uLen = length(u)

  if (uLen < 1e-12) {
    u = pB // Try pB if pA was (anti)parallel to n
    u = sub(u, mul(nNorm, dot(u, nNorm)))
    uLen = length(u)
    if (uLen < 1e-12) {
      return null // pA and pB are both (anti)parallel to n? Should not happen if nLen > 0.
    }
  }
  u = mul(u, 1 / uLen) // Normalize u

  const v = cross(nNorm, u) // v is perpendicular to both n and u

  // Point on great circle: p(θ) = cos(θ)*u + sin(θ)*v
  // We need z-component: p(θ).z = sin(latC)
  // cos(θ)*u.z + sin(θ)*v.z = sin(latC)
  const uz = u[2]
  const vz = v[2]
  const target = Math.sin(toRad(latC))

  // Solve A*cos(θ) + B*sin(θ) = C  (where A=uz, B=vz, C=target)
  const R = Math.sqrt(uz * uz + vz * vz)

  if (R < 1e-12) {
    // Great circle is on the equator (uz=0, vz=0)
    if (Math.abs(target) < 1e-12) { // latC is 0
      // All points on equator intersect. Return two antipodal points.
      const p0 = u
      const p180 = mul(u, -1)
      return [normalizeLon(cartesianToLonLat(p0)[0]), normalizeLon(cartesianToLonLat(p180)[0])]
    }
    return null // Equator circle never intersects non-zero latitude
  }

  // Check if solution exists
  const ratio = target / R
  if (Math.abs(ratio) > 1 + 1e-12) {
    return null // Target latitude is "higher" than the circle's max latitude
  }

  const ratioClamp = Math.max(-1, Math.min(1, ratio))
  const phi = Math.atan2(vz, uz) // Phase angle

  // θ - φ = ±acos(C/R) => θ = φ ± acos(C/R)
  const delta = Math.acos(ratioClamp)
  const theta1 = phi + delta
  const theta2 = phi - delta

  // Get intersection points in cartesian
  const p1_cart = add(mul(u, Math.cos(theta1)), mul(v, Math.sin(theta1)))
  const p2_cart = add(mul(u, Math.cos(theta2)), mul(v, Math.sin(theta2)))

  // Convert to lon/lat
  const lon1 = normalizeLon(cartesianToLonLat(p1_cart)[0])
  const lon2 = normalizeLon(cartesianToLonLat(p2_cart)[0])

  return [lon1, lon2]
}

// Input: p1LonLat, p2LonLat as [lon, lat]
// Output: Array of intersection points [lon, lat]
export function getIntersections(p1LonLat, p2LonLat, GAME_BOUNDS = null) {
  const bisectionLine = getBisectionLine(p1LonLat, p2LonLat)
  // Reads from module-level variable
  let northBoundIntersections = getIntersectionWithConstantLat(bisectionLine[0], bisectionLine[1], CURRENT_GAME_BOUNDS.north)
  let southBoundIntersections = getIntersectionWithConstantLat(bisectionLine[0], bisectionLine[1], CURRENT_GAME_BOUNDS.south)

  if (GAME_BOUNDS !== null) {
    CURRENT_GAME_BOUNDS = GAME_BOUNDS;
    CURRENT_GAME_BOUNDS_CORNERS_LONLAT = [
      [GAME_BOUNDS.west, GAME_BOUNDS.south], // SW
      [GAME_BOUNDS.west, GAME_BOUNDS.north], // NW
      [GAME_BOUNDS.east, GAME_BOUNDS.north], // NE
      [GAME_BOUNDS.east, GAME_BOUNDS.south]  // SE
    ];
  }
  // These are intersections with the *great circles* forming the boundary
  const rectangleIntersections = getBoundaryIntersections(p1LonLat, p2LonLat);

  let intersections = [];

  // Add intersections with North latitude line
  // Reads from module-level variable
  if (!northBoundIntersections) northBoundIntersections = []
  for (const lon of northBoundIntersections) {
    if (lon >= CURRENT_GAME_BOUNDS.west && lon <= CURRENT_GAME_BOUNDS.east) {
      intersections.push([lon, CURRENT_GAME_BOUNDS.north])
    }
  }

  // Add intersections with South latitude line
  // Reads from module-level variable
  if (!southBoundIntersections) southBoundIntersections = []
  for (const lon of southBoundIntersections) {
    if (lon >= CURRENT_GAME_BOUNDS.west && lon <= CURRENT_GAME_BOUNDS.east) {
      intersections.push([lon, CURRENT_GAME_BOUNDS.south])
    }
  }

  const e = 1e-6;

  // Add intersections with West/East longitude lines
  // Reads from module-level variable
  for (const inter of rectangleIntersections) {
    if (Math.abs(inter[0] - CURRENT_GAME_BOUNDS.west) < e ||
      Math.abs(inter[0] - CURRENT_GAME_BOUNDS.east) < e) {
      if (inter[1] >= CURRENT_GAME_BOUNDS.south && inter[1] <= CURRENT_GAME_BOUNDS.north) {
        intersections.push(inter)
      }
    }
  }

  return intersections
}

// Input: p1LonLat, p2LonLat as [lon, lat]
// Output: Array of 0s (p1 closer) or 1s (p2 closer) for each corner
function checkSideOfCorner(p1LonLat, p2LonLat) {
  let sides = [];

  // Convert input points to cartesian once
  const p1Cartesian = lonLatToCartesian(p1LonLat[0], p1LonLat[1]);
  const p2Cartesian = lonLatToCartesian(p2LonLat[0], p2LonLat[1]);

  // Iterate over [lon, lat] corners
  // Reads from module-level variable
  for (const cornerLonLat of CURRENT_GAME_BOUNDS_CORNERS_LONLAT) {
    const cornerCartesian = lonLatToCartesian(cornerLonLat[0], cornerLonLat[1]);

    // Compare distances in cartesian space
    const dist1 = sphericalDistance(p1Cartesian, cornerCartesian);
    const dist2 = sphericalDistance(p2Cartesian, cornerCartesian);

    if (dist1 < dist2) {
      sides.push(0); // p1 is closer
    } else if (dist2 < dist1) {
      sides.push(1); // p2 is closer
    }
    // Note: if distances are equal, no side is pushed. This matches original logic.
  }
  return sides;
}

/**
 * Builds a polygon for one side of the bisection line.
 * @param {object} GAME_BOUNDS - An object with {north, south, east, west} properties.
 * @param {Array<Array<number>>} pts - An array of two [lon, lat] points: [p1LonLat, p2LonLat].
 * @param {boolean} side - The side to build the polygon for (false/0 for p1, true/1 for p2).
 * @returns {Array<Array<number>>} An array of [lon, lat] coordinates forming the polygon.
 */
export function buildPolygon(GAME_BOUNDS, pts, side, resolution=100) {
  // 1. Set the module-level (private) variables
  CURRENT_GAME_BOUNDS = GAME_BOUNDS;
  CURRENT_GAME_BOUNDS_CORNERS_LONLAT = [
    [CURRENT_GAME_BOUNDS.west, CURRENT_GAME_BOUNDS.south], // SW
    [CURRENT_GAME_BOUNDS.west, CURRENT_GAME_BOUNDS.north], // NW
    [CURRENT_GAME_BOUNDS.east, CURRENT_GAME_BOUNDS.north], // NE
    [CURRENT_GAME_BOUNDS.east, CURRENT_GAME_BOUNDS.south]  // SE
  ];

  const p1LonLat = pts[0];
  const p2LonLat = pts[1];

  // 2. Call internal functions that use the module-level variables
  const boundaryIntersections = getIntersections(p1LonLat, p2LonLat);
  
  // Guard against cases where the bisection line doesn't produce two intersections
  if (boundaryIntersections.length < 2) {
      console.error("buildPolygon: Could not find two boundary intersections. Returning empty polygon.");
      return [];
  }

  const cornerSides = checkSideOfCorner(p1LonLat, p2LonLat);
  const bisectionCoords = getLineCoordsBetween(boundaryIntersections[0], boundaryIntersections[1], resolution);

  let cornerCoords = [];

  for (let i in cornerSides) {
    // `side` (boolean) will be coerced to 0 (false) or 1 (true) by the `==` operator
    if (cornerSides[i] != side) {
      cornerCoords.push([CURRENT_GAME_BOUNDS_CORNERS_LONLAT[i]]);
    }
  }

  let polygonCoords = [];
  polygonCoords = polygonCoords.concat(bisectionCoords);

  let added = true;
  while (added) {
    added = false;
    for (let j = cornerCoords.length - 1; j >= 0; j--) {
      const corner = cornerCoords[j];
      const last = polygonCoords[polygonCoords.length - 1];
      if (
        Math.abs(corner[0][0] - last[0]) < 1e-6 ||
        Math.abs(corner[0][1] - last[1]) < 1e-6
      ) {
        polygonCoords.push(corner[0]);
        cornerCoords.splice(j, 1);
        added = true;
      }
    }
  }

  return polygonCoords;
}