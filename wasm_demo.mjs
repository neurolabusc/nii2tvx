// node wasm_demo.mjs lesion.nii[.gz] atlas.tvx   (build first: make wasm)
// Prints the same TSV as the native tool. Shows the whole WASM surface.
import { readFileSync } from "node:fs";
import { gunzipSync } from "node:zlib";
import { basename } from "node:path";
import createModule from "./nii2tvx.mjs";

const [lesionPath, tvxPath] = process.argv.slice(2);
const M = await createModule();

function toHeap(bytes) { // copy a Node buffer into WASM memory; caller frees
	const ptr = M._malloc(bytes.length);
	M.HEAPU8.set(bytes, ptr);
	return ptr;
}

const tvxBytes = readFileSync(tvxPath);
const tvx = M._tvx_open(toHeap(tvxBytes), tvxBytes.length); // tvx_open owns the buffer
if (!tvx) process.exit(1);
const ntract = M._tvx_ntract(tvx);
const names = Array.from({ length: ntract }, (_, k) => M.UTF8ToString(M._tvx_name(tvx, k)));

let nii = readFileSync(lesionPath);
if (nii[0] === 0x1f && nii[1] === 0x8b) nii = gunzipSync(nii); // mask_open wants uncompressed NIfTI
const niiPtr = toHeap(nii);
const mask = M._mask_open(niiPtr, nii.length);
M._free(niiPtr);
if (!mask) process.exit(1);

const fracs = [];
for (let k = 0; k < ntract; k++) {
	const f = M._tvx_query(tvx, k, mask);
	if (f < 0) process.exit(1); // grid mismatch or corrupt file; reason was printed to stderr
	fracs.push(f);
}
M._mask_close(mask);
M._tvx_close(tvx);

console.log(["id", ...names].join("\t"));
console.log([basename(lesionPath).replace(/\.nii(\.gz)?$/, ""), ...fracs.map((f) => Number(f.toPrecision(6)))].join("\t"));
