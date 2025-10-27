let GAME_BOUNDS = [];

function getBisectorLinePoints(p1_xy, p2_xy) {
    const rect = { x: GAME_BOUNDS.west, y: GAME_BOUNDS.south, width: GAME_BOUNDS.east - GAME_BOUNDS.west, height: GAME_BOUNDS.north - GAME_BOUNDS.south };
    const dx = p2_xy.x - p1_xy.x;
    const dy = p2_xy.y - p1_xy.y;
    if (Math.abs(dx) < 1e-9 && Math.abs(dy) < 1e-9) return [];
    const midPoint = { x: (p1_xy.x + p2_xy.x) / 2, y: (p1_xy.y + p2_xy.y) / 2 };
    const a = dx; const b = dy; const c = -dx * midPoint.x - dy * midPoint.y;
    const intersections = [];
    if (Math.abs(a) > 1e-9) { const x_top = (-c - b * (rect.y + rect.height)) / a; if (x_top >= rect.x && x_top <= rect.x + rect.width) intersections.push({ y: rect.y + rect.height, x: x_top }); }
    if (Math.abs(a) > 1e-9) { const x_bottom = (-c - b * rect.y) / a; if (x_bottom >= rect.x && x_bottom <= rect.x + rect.width) intersections.push({ y: rect.y, x: x_bottom }); }
    if (Math.abs(b) > 1e-9) { const y_left = (-c - a * rect.x) / b; if (y_left >= rect.y && y_left <= rect.y + rect.height) intersections.push({ y: y_left, x: rect.x }); }
    if (Math.abs(b) > 1e-9) { const y_right = (-c - a * (rect.x + rect.width)) / b; if (y_right >= rect.y && y_right <= rect.y + rect.height) intersections.push({ y: y_right, x: rect.x + rect.width }); }
    const uniquePoints = []; const seen = new Set();
    for (const p of intersections) { const key = `${p.x.toFixed(7)},${p.y.toFixed(7)}`; if (!seen.has(key)) { uniquePoints.push(p); seen.add(key); } }
    return uniquePoints.length >= 2 ? [[uniquePoints[0].y, uniquePoints[0].x], [uniquePoints[1].y, uniquePoints[1].x]] : [];
}

function getRectangleSlice(rect, pointA, pointB) {
    const dx = pointB.x - pointA.x; const dy = pointB.y - pointA.y;
    if (Math.abs(dx) < 1e-9 && Math.abs(dy) < 1e-9) { return [{ x: rect.x, y: rect.y }, { x: rect.x + rect.width, y: rect.y }, { x: rect.x + rect.width, y: rect.y + rect.height }, { x: rect.x, y: rect.y + rect.height }]; }
    const midPoint = { x: (pointA.x + pointB.x) / 2, y: (pointA.y + pointB.y) / 2 };
    const line = { a: dx, b: dy, c: -dx * midPoint.x - dy * midPoint.y };
    const sideOfA = line.a * pointA.x + line.b * pointA.y + line.c;
    const isInside = (p) => (line.a * p.x + line.b * p.y + line.c) * sideOfA >= 0;
    const subjectPolygon = [{ x: rect.x, y: rect.y }, { x: rect.x + rect.width, y: rect.y }, { x: rect.x + rect.width, y: rect.y + rect.height }, { x: rect.x, y: rect.y + rect.height }];
    const outputPolygon = []; let S = subjectPolygon[subjectPolygon.length - 1];
    for (const E of subjectPolygon) {
        const sIsInside = isInside(S); const eIsInside = isInside(E);
        if (eIsInside) { if (!sIsInside) outputPolygon.push(getIntersection(S, E, line)); outputPolygon.push(E); } else if (sIsInside) { outputPolygon.push(getIntersection(S, E, line)); }
        S = E;
    } return outputPolygon;
}

function getIntersection(p1, p2, line) {
    const { a, b, c } = line; const dx = p2.x - p1.x; const dy = p2.y - p1.y;
    const denominator = a * dx + b * dy; if (Math.abs(denominator) < 1e-9) return p1;
    const t = -(a * p1.x + b * p1.y + c) / denominator; return { x: p1.x + t * dx, y: p1.y + t * dy };
}

function isLeft(x1, y1, x2, y2, px, py) {
    return (x2 - x1) * (py - y1) - (px - x1) * (y2 - y1)
}

/**
 * Builds a polygon for one side of the bisection line.
 * @param {object} GAME_BOUNDS - An object with {north, south, east, west} properties.
 * @param {Array<Array<number>>} pts - An array of two [lon, lat] points: [p1LonLat, p2LonLat].
 * @param {boolean} side - The side to build the polygon for (false/0 for p1, true/1 for p2).
 * @returns {Array<Array<number>>} An array of [lon, lat] coordinates forming the polygon.
 */
export function createEliminationZone(prevBar, currentBar, GAME_BOUNDS, isWarmer) {
    GAME_BOUNDS = GAME_BOUNDS;
    
    const rect = { x: GAME_BOUNDS.west, y: GAME_BOUNDS.south, width: GAME_BOUNDS.east - GAME_BOUNDS.west, height: GAME_BOUNDS.north - GAME_BOUNDS.south };
    const pointToEliminate = isWarmer ? prevBar : currentBar; const pointToKeep = isWarmer ? currentBar : prevBar;
    const pointA_xy = { x: pointToEliminate.lng, y: pointToEliminate.lat }; const pointB_xy = { x: pointToKeep.lng, y: pointToKeep.lat };
    const polygonVertices = getRectangleSlice(rect, pointA_xy, pointB_xy); if (polygonVertices.length < 3) return null;
    const polygon = polygonVertices.map(obj => [obj.x, obj.y]);
    return polygon;
}

// HACKER: Find if point is inside polygon (elimination zone)
export function isPointInPolygon(point, polygon) {
    let wn = 0
    const px = point.lng
    const py = point.lat

    for (let i = 0; i < polygon.length; i++) {
        const { x: x1, y: y1 } = polygon[i]
        const { x: x2, y: y2 } = polygon[(i + 1) % polygon.length]

        if (y1 <= py) {
            if (y2 > py && isLeft(x1, y1, x2, y2, px, py) > 0) wn++
        } else {
            if (y2 <= py && isLeft(x1, y1, x2, y2, px, py) < 0) wn--
        }
    }

    return wn !== 0
}

