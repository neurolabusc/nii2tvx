// gcc -O3 nii2tvx.c -o nii2tvx -lz -lm
// gcc -O1 -g -fsanitize=address -fno-omit-frame-pointer nii2tvx.c -o nii2tvx -lz -lm
#include <stdio.h>
//#include <ctype.h>
#include <stdlib.h>
//#include <stdint.h>
#include <math.h>
#include <string.h>
#include "nifti1.h"
#include <zlib.h>
#include<stdbool.h>
#ifdef _MSC_VER
	#define F_OK 0 
#endif

#define STR_HELPER(x) #x
#define STR(x) STR_HELPER(x)

#if defined(__ICC) || defined(__INTEL_COMPILER)
	#define kCCsuf " IntelCC" STR(__INTEL_COMPILER)
#elif defined(_MSC_VER)
	#define kCCsuf " MSC" STR(_MSC_VER)
#elif defined(__clang__)
	#define kCCsuf " Clang" STR(__clang_major__) "." STR(__clang_minor__) "." STR(__clang_patchlevel__)
#elif defined(__GNUC__) || defined(__GNUG__)
	#define kCCsuf " GCC" STR(__GNUC__) "." STR(__GNUC_MINOR__) "." STR(__GNUC_PATCHLEVEL__)
#else
	#define kCCsuf " CompilerNA"
#endif
#if defined(__arm__) || defined(__ARM_ARCH)
	#define kCPUsuf " ARM"
#elif defined(__x86_64)
	#define kCPUsuf " x86-64"
#else
	#define kCPUsuf " " //unknown CPU
#endif
#define kdate "v1.0.20240101"

bool is_ext(const char *filename, const char *extension) {
	size_t filename_len = strlen(filename);
	size_t extension_len = strlen(extension);
	// Check if the filename is long enough to contain the extension
	if (filename_len >= extension_len) {
		const char *file_extension = filename + (filename_len - extension_len);
		// Compare the file extension with the specified extension
		if (strcmp(file_extension, extension) == 0) {
			return true;
		}
	}
	return false;
}

void strip_ext(char *fname){
	char *end = fname + strlen(fname);
	while (end > fname && *end != '.' && *end != '\\' && *end != '/') {
		--end;
	}
	if ((end > fname && *end == '.') &&
		(*(end - 1) != '\\' && *(end - 1) != '/')) {
		*end = '\0';
	}
}

void strip_ext2(char *fname) {
	char *end = fname + strlen(fname);
	// Check if the file has ".gz" extension
	if (end >= fname + 3 && strcmp(end - 3, ".gz") == 0)
		end -= 4; // Move end back by 4 characters to exclude ".gz"
	while (end > fname && *end != '.' && *end != '\\' && *end != '/') {
		--end;
	}
	if ((end > fname && *end == '.') &&
		(*(end - 1) != '\\' && *(end - 1) != '/')) {
		*end = '\0';
	}
}

const char *basenamex(const char *path) {
	const char *s = strrchr(path, '/');
	return s ? s + 1 : path;
}

float * load_nii(const char *fnm, nifti_1_header * hdr) {
	char imgnm[768], hdrnm[768], basenm[768], ext[768] = "";
	strcpy(basenm, fnm);
	strcpy(imgnm, fnm);
	strcpy(hdrnm, fnm);
	strip_ext(basenm); // ~/file.nii -> ~/file
	if (strlen(fnm) > strlen(basenm))
		strcpy(ext, fnm + strlen(basenm));
	if (strstr(ext, ".hdr")) {
		strcpy(imgnm, basenm);
		strcat(imgnm, ".img");
	}
	if (strstr(ext, ".img")) {
		strcpy(hdrnm, basenm);
		strcat(hdrnm, ".hdr");
	}
	if( access( hdrnm, F_OK ) != 0 ) {
		printf("Unable to find a file named %s\n", hdrnm);
		return NULL;
	}
	if( access( imgnm, F_OK ) != 0 ) {
		printf("Unable to find a file named %s\n", imgnm);
		return NULL;
	}
	bool isGz = false;
	if (strstr(ext, ".gz")) {
		gzFile fgz = gzopen(hdrnm, "r");
		if (! fgz) {
			printf("gzopen error %s\n", hdrnm);
			return NULL;
		}
		int bytes_read = gzread(fgz, hdr, sizeof(nifti_1_header));
		gzclose(fgz);
		isGz = true;
		if (bytes_read < sizeof(nifti_1_header)) return NULL;
	}
	if (!isGz) {
		FILE *fp = fopen(hdrnm,"rb");
		if (fp == NULL) {
			printf("Unable to open %s\n", hdrnm);
			return NULL;
		}
		size_t bytes_read = fread(hdr, sizeof(nifti_1_header), 1, fp);
		fclose(fp);
		if (bytes_read <= 0) {
			printf("Unable to read %s\n", hdrnm);
			return NULL;
		}
	}
	uint16_t sig = 348;
	uint16_t fwd = hdr->sizeof_hdr;
	uint16_t rev = (fwd & 0xff) << 8 | ((fwd & 0xff00) >> 8);
	if (rev == sig) {
		printf("Demo only reads native endian NIfTI (solution: use niimath)\n");
		return NULL;
	}
	if (fwd != sig) {
		printf("Only compiled to read uncompressed NIfTI (solution: 'gzip -d \"%s\"')\n", hdrnm);
		return NULL;
	}
	if ((hdr->datatype != DT_UINT8) && (hdr->datatype != DT_UINT16) && (hdr->datatype != DT_INT16) && (hdr->datatype != DT_FLOAT32)) {
		printf("Demo does not support this data type (solution: use niimath)\n");
		return NULL;
	}
	int nvox = hdr->dim[1]*hdr->dim[2]*hdr->dim[3];
	if (hdr->scl_slope == 0.0) hdr->scl_slope = 1.0;
	float * img32 = (float *) malloc(nvox*sizeof(float));
	int bpp = 2;
	if (hdr->datatype == DT_UINT8)
		bpp = 1;
	if (hdr->datatype == DT_FLOAT32)
		bpp = 4;
	void * imgRaw = (void *) malloc(nvox*bpp);
	if (isGz) {
		gzFile fgz = gzopen(hdrnm, "r");
		if (! fgz)
			return NULL;
		if (hdr->vox_offset > 0) { //skip header
			int nskip = round(hdr->vox_offset);
			void * skip = (void *) malloc(nskip);
			int bytes_read = gzread (fgz, skip, nskip);
			free(skip);
			if (bytes_read != nskip)
				return NULL;
		}
		int bytes_read = gzread(fgz, imgRaw, nvox*bpp);
		if (bytes_read != (nvox*bpp))
			return NULL;
		gzclose(fgz);
	} else {
		FILE *fp = fopen(imgnm,"rb");
		if (fp == NULL)
			return NULL;
		fseek(fp, (int)hdr->vox_offset, SEEK_SET);
		size_t sz = fread(imgRaw, nvox*bpp, 1, fp);
		fclose(fp);
		if (sz <= 0) return NULL;
	}
	if (hdr->datatype == DT_UINT8) {
		uint8_t * img8 = (uint8_t *) imgRaw;
		for (int i = 0; i < nvox; i++)
			img32[i] = (img8[i] * hdr->scl_slope) + hdr->scl_inter;
	} else if (hdr->datatype == DT_UINT16) {
		uint16_t * img16 = (uint16_t *) imgRaw;
		for (int i = 0; i < nvox; i++)
			img32[i] = (img16[i] * hdr->scl_slope) + hdr->scl_inter;
	} else if (hdr->datatype == DT_INT16) {
		int16_t * img16 = (int16_t *) imgRaw;
		for (int i = 0; i < nvox; i++)
			img32[i] = (img16[i] * hdr->scl_slope) + hdr->scl_inter;
	} else {
		float * img32w = (float *) imgRaw;
		for (int i = 0; i < nvox; i++)
			img32[i] = (img32w[i] * hdr->scl_slope) + hdr->scl_inter;
	}
	free(imgRaw);
	return img32;
}

uint8_t * load_nii_mask(const char *fnm, nifti_1_header * hdr) {
	//image data is binarized as 0 or 1 (not zero)
	float * img32 = load_nii(fnm, hdr);
	if (img32 == NULL) {
		return NULL;
	}
	int nvox = hdr->dim[1]*hdr->dim[2]*hdr->dim[3];
	uint8_t * img = (uint8_t*)calloc(nvox, sizeof(uint8_t));
	
	for (int i = 0; i < nvox; i++)
		if (img32[i] != 0.0)
			img[i] = 1;
	free(img32);
	return img;
}

typedef struct {/** 4x4 matrix struct **/
	float m[4][4];
} mat44;

typedef struct {
	int32_t v[3];
} vec3i;

typedef struct { /** x4 vector struct **/
	float v[4];
} vec4;

mat44 nifti_mat44_mul( mat44 A , mat44 B ) {
	mat44 C ; int i,j,k;
	for( i=0 ; i < 4 ; i++ )
		for( j=0 ; j < 4 ; j++ ) {
			C.m[i][j] = 0.0;
			for( k=0; k < 4; k++ )
				C.m[i][j] += A.m[i][k] * B.m[k][j];
		}
	return C;
}

vec3i nifti_xyz_mat44_mul(float x, float y, float z, mat44 m) { //multiply vector * 4x4matrix
	vec4 v = {{x, y, z, 1}};
	vec4 vO;
	for (int i = 0; i < 4; i++) { //multiply Pcrs * m
		vO.v[i] = 0;
		for (int j = 0; j < 4; j++)
			vO.v[i] += m.m[i][j] * v.v[j];
	}
	
	vec3i vOi = {{round(vO.v[0]), round(vO.v[1]), round(vO.v[2])}};
	return vOi;
}

mat44 nifti_mat44_inverse( mat44 R ) {
	double r11,r12,r13,r21,r22,r23,r31,r32,r33,v1,v2,v3 , deti;
	mat44 Q;
	/* INPUT MATRIX IS: */
	r11 = R.m[0][0]; r12 = R.m[0][1]; r13 = R.m[0][2]; /* [ r11 r12 r13 v1 ] */
	r21 = R.m[1][0]; r22 = R.m[1][1]; r23 = R.m[1][2]; /* [ r21 r22 r23 v2 ] */
	r31 = R.m[2][0]; r32 = R.m[2][1]; r33 = R.m[2][2]; /* [ r31 r32 r33 v3 ] */
	v1 = R.m[0][3]; v2 = R.m[1][3]; v3 = R.m[2][3]; /* [ 0 0 0 1 ] */
	deti = r11*r22*r33-r11*r32*r23-r21*r12*r33
		 +r21*r32*r13+r31*r12*r23-r31*r22*r13;
	if( deti != 0.0l ) deti = 1.0l / deti;
	Q.m[0][0] = (float)( deti*( r22*r33-r32*r23) );
	Q.m[0][1] = (float)( deti*(-r12*r33+r32*r13) );
	Q.m[0][2] = (float)( deti*( r12*r23-r22*r13) );
	Q.m[0][3] = (float)( deti*(-r12*r23*v3+r12*v2*r33+r22*r13*v3
					 -r22*v1*r33-r32*r13*v2+r32*v1*r23) );
	Q.m[1][0] = (float)( deti*(-r21*r33+r31*r23) );
	Q.m[1][1] = (float)( deti*( r11*r33-r31*r13) );
	Q.m[1][2] = (float)( deti*(-r11*r23+r21*r13) );
	Q.m[1][3] = (float)( deti*( r11*r23*v3-r11*v2*r33-r21*r13*v3
					 +r21*v1*r33+r31*r13*v2-r31*v1*r23) );
	Q.m[2][0] = (float)( deti*( r21*r32-r31*r22) );
	Q.m[2][1] = (float)( deti*(-r11*r32+r31*r12) );
	Q.m[2][2] = (float)( deti*( r11*r22-r21*r12) );
	Q.m[2][3] = (float)( deti*(-r11*r22*v3+r11*r32*v2+r21*r12*v3
					 -r21*r32*v1-r31*r12*v2+r31*r22*v1) );
	Q.m[3][0] = Q.m[3][1] = Q.m[3][2] = 0.0l;
	Q.m[3][3] = (deti == 0.0l) ? 0.0l : 1.0l ; /* failure flag if deti == 0 */
	return Q;
}

mat44 make_mat44(
	float m00, float m01, float m02, float m03,
	float m10, float m11, float m12, float m13,
	float m20, float m21, float m22, float m23)
{
	mat44 m;
	m.m[0][0] = m00;
	m.m[0][1] = m01;
	m.m[0][2] = m02;
	m.m[0][3] = m03;
	m.m[1][0] = m10;
	m.m[1][1] = m11;
	m.m[1][2] = m12;
	m.m[1][3] = m13;
	m.m[2][0] = m20;
	m.m[2][1] = m21;
	m.m[2][2] = m22;
	m.m[2][3] = m23;
	m.m[3][0] = 0;
	m.m[3][1] = 0;
	m.m[3][2] = 0;
	m.m[3][3] = 1;
	return m;
}

mat44 sform(nifti_1_header * hdr) {
	return make_mat44(
		hdr->srow_x[0], hdr->srow_x[1], hdr->srow_x[2], hdr->srow_x[3],
		hdr->srow_y[0], hdr->srow_y[1], hdr->srow_y[2], hdr->srow_y[3],
		hdr->srow_z[0], hdr->srow_z[1], hdr->srow_z[2], hdr->srow_z[3]);
}

#pragma pack(2)
struct trk_header { //always little endian
	char id_string[6];
	int16_t dim[3];
	float voxel_size[3];
	float origin[3];
	int16_t n_scalars;
	char scalar_name[10][20];
	int16_t n_properties;
	char property_name[10][20];
	mat44 vox_to_ras;
	char reserved[444];
	char voxel_order[4];
	char pad2[4];
	float image_orientation_patient[6];
	char pad1[2];
	uint8_t invert_x, invert_y,swap_xy, swap_yz, swap_zx;
	int32_t n_count, version, hdr_size;
};
#pragma pack()
typedef struct trk_header trk_header;

#pragma pack(2)
struct tvx_header { //always little endian
	uint32_t signature; //must be 175666804 "tvx\n"
	uint32_t dim[3]; //correspond to NIfTI dim[1..3]
	float srow_x[4] ; // 1st row affine transform.
	float srow_y[4] ; // 2nd row affine transform.
	float srow_z[4] ; // 3rd row affine transform.
	uint32_t noffset; //number of offsets (streamlines +1)
	uint32_t nvoxel; //number of voxels stored for all streamlines
};
#pragma pack()
typedef struct tvx_header tvx_header;


#define kSigRaw 0x0A787674   // "tvx\n": uint32 voxel indices
#define kSigDelta 0x0A647674 // "tvd\n": one byte per voxel step, see tvx_format.md

typedef struct {
	tvx_header h;
	uint32_t *offsets, *verts;
} tvx_t;

// ---- delta codec: 26 neighbour codes (0..26, 13 unused) + escape 27 followed by zigzag LEB128 delta
static void neighbour_table(int64_t tab[27], int nx, int ny) {
	for (int c = 0; c < 27; c++)
		tab[c] = (c % 3 - 1) + (int64_t)(c / 3 % 3 - 1) * nx + (int64_t)(c / 9 - 1) * nx * ny;
}

static size_t delta_encode(const tvx_t *t, uint8_t *out) {
	int64_t tab[27];
	neighbour_table(tab, t->h.dim[0], t->h.dim[1]);
	uint8_t *p = out;
	int64_t prev = 0;
	for (uint32_t i = 0; i < t->h.nvoxel; i++) {
		int64_t d = (int64_t)t->verts[i] - prev;
		prev = t->verts[i];
		int c = 0;
		while (c < 27 && (c == 13 || tab[c] != d))
			c++;
		*p++ = c;
		if (c == 27) {
			uint64_t z = ((uint64_t)d << 1) ^ (uint64_t)(d >> 63);
			for (; z >= 0x80; z >>= 7)
				*p++ = 0x80 | (z & 0x7F);
			*p++ = z;
		}
	}
	return p - out;
}

static void delta_decode(tvx_t *t, const uint8_t *in) {
	int64_t tab[27];
	neighbour_table(tab, t->h.dim[0], t->h.dim[1]);
	int64_t prev = 0;
	for (uint32_t i = 0; i < t->h.nvoxel; i++) {
		int c = *in++;
		if (c < 27)
			prev += tab[c];
		else {
			uint64_t z = 0;
			int shift = 0;
			for (;; shift += 7) {
				uint8_t b = *in++;
				z |= (uint64_t)(b & 0x7F) << shift;
				if (!(b & 0x80)) break;
			}
			prev += (int64_t)(z >> 1) ^ -(int64_t)(z & 1);
		}
		t->verts[i] = prev;
	}
}

// ---- conversion: streamline vertices (mm) -> voxel index runs
typedef struct {
	tvx_t t;
	size_t vcap, ocap, nvert; // nvert counts input vertices
	int prev;
	bool has_prev;
	float pv[3]; // previous vertex, continuous voxel coords
	mat44 inv;   // mm (or TRK voxmm) -> voxel
} tvx_writer;

static void writer_init(tvx_writer *w, nifti_1_header *hdr, mat44 inv) {
	memset(w, 0, sizeof(*w));
	w->t.h.signature = kSigRaw;
	for (int i = 0; i < 3; i++) w->t.h.dim[i] = hdr->dim[i + 1];
	for (int i = 0; i < 4; i++) {
		w->t.h.srow_x[i] = hdr->srow_x[i];
		w->t.h.srow_y[i] = hdr->srow_y[i];
		w->t.h.srow_z[i] = hdr->srow_z[i];
	}
	w->inv = inv;
	w->prev = -1;
	w->ocap = 1 << 16;
	w->t.offsets = malloc(w->ocap * sizeof(uint32_t));
	w->t.offsets[0] = 0;
	w->t.h.noffset = 1;
}

static void emit(tvx_writer *w, int x, int y, int z) {
	uint32_t *dim = w->t.h.dim;
	if (x < 0 || y < 0 || z < 0 || x >= (int)dim[0] || y >= (int)dim[1] || z >= (int)dim[2])
		return;
	int vxl = x + y * dim[0] + z * dim[0] * dim[1];
	if (vxl == w->prev)
		return;
	w->prev = vxl;
	if (w->t.h.nvoxel == w->vcap) {
		w->vcap = w->vcap ? 2 * w->vcap : 1 << 20;
		w->t.verts = realloc(w->t.verts, w->vcap * sizeof(uint32_t));
	}
	w->t.verts[w->t.h.nvoxel++] = vxl;
}

// Amanatides & Woo grid traversal: every voxel the segment a->b crosses.
// Voxel v spans [v-0.5, v+0.5), matching nearest-voxel rounding.
static void raster(tvx_writer *w, const float *a, const float *b) {
	int v[3], step[3];
	float tMax[3], tDelta[3];
	for (int k = 0; k < 3; k++) {
		v[k] = (int)floorf(a[k] + 0.5f);
		float d = b[k] - a[k];
		step[k] = (d > 0) - (d < 0);
		tDelta[k] = step[k] ? fabsf(1.0f / d) : INFINITY;
		tMax[k] = step[k] ? (v[k] + 0.5f * step[k] - a[k]) / d : INFINITY;
	}
	emit(w, v[0], v[1], v[2]);
	for (;;) {
		int k = tMax[0] < tMax[1] ? (tMax[0] < tMax[2] ? 0 : 2) : (tMax[1] < tMax[2] ? 1 : 2);
		if (tMax[k] > 1.0f)
			break;
		v[k] += step[k];
		tMax[k] += tDelta[k];
		emit(w, v[0], v[1], v[2]);
	}
}

static void add_vertex(tvx_writer *w, const float *xyz) {
	float p[3];
	for (int i = 0; i < 3; i++)
		p[i] = xyz[0] * w->inv.m[i][0] + xyz[1] * w->inv.m[i][1] + xyz[2] * w->inv.m[i][2] + w->inv.m[i][3];
	w->nvert++;
	raster(w, w->has_prev ? w->pv : p, p);
	memcpy(w->pv, p, sizeof(p));
	w->has_prev = true;
}

static void end_streamline(tvx_writer *w) {
	uint32_t n = w->t.h.noffset;
	if (w->t.h.nvoxel > w->t.offsets[n - 1]) { // drop streamlines with no in-volume voxel
		if (n == w->ocap)
			w->t.offsets = realloc(w->t.offsets, (w->ocap *= 2) * sizeof(uint32_t));
		w->t.offsets[n] = w->t.h.nvoxel;
		w->t.h.noffset++;
	}
	w->has_prev = false;
	w->prev = -1;
}

static void write_tvx(const char *fnm, tvx_writer *w, bool delta) {
	char outnm[768];
	strcpy(outnm, fnm);
	strip_ext(outnm);
	strcat(outnm, ".tvx");
	FILE *fp = fopen(outnm, "wb");
	if (!fp) {
		printf("Unable to write %s\n", outnm);
		exit(EXIT_FAILURE);
	}
	tvx_t *t = &w->t;
	t->h.signature = delta ? kSigDelta : kSigRaw;
	fwrite(&t->h, sizeof(t->h), 1, fp);
	fwrite(t->offsets, sizeof(uint32_t), t->h.noffset, fp);
	if (delta) {
		uint8_t *buf = malloc((size_t)t->h.nvoxel * 11); // worst case: escape + 10-byte varint
		fwrite(buf, 1, delta_encode(t, buf), fp);
		free(buf);
	} else
		fwrite(t->verts, sizeof(uint32_t), t->h.nvoxel, fp);
	fclose(fp);
	printf("%s\t%u\tstreamlines\t%zu\tvertices\t%u\tvoxels\n", outnm, t->h.noffset - 1, w->nvert, t->h.nvoxel);
	free(t->offsets);
	free(t->verts);
}

static int load_trk(const char *fnm, nifti_1_header *hdr, bool delta) {
	//https://trackvis.org/docs/?subsect=fileformat
	FILE *fp = fopen(fnm, "rb");
	if (fp == NULL)
		return EXIT_FAILURE;
	trk_header thdr;
	if (fread(&thdr, sizeof(trk_header), 1, fp) != 1 || thdr.hdr_size != sizeof(trk_header) || thdr.version != 2 || thdr.n_count == 0) {
		printf("Unable to read TRK header %d %d\n", thdr.hdr_size, thdr.version);
		if (thdr.n_count == 0)
			printf("unable to read TRK with implicit n_count (hint: convert to TCK with tff_convert_tractogram.py)\n");
		fclose(fp);
		return EXIT_FAILURE;
	}
	if (thdr.vox_to_ras.m[3][3] == 0.0) {
		printf("TRK vox_to_ras not set\n");
		fclose(fp);
		return EXIT_FAILURE;
	}
	// TRK vertices are in voxmm with corner origin: scale to voxels, shift half a voxel, then vox_to_ras
	mat44 zoomMat = make_mat44(
		1.0 / thdr.voxel_size[0], 0, 0, -0.5,
		0, 1.0 / thdr.voxel_size[1], 0, -0.5,
		0, 0, 1.0 / thdr.voxel_size[2], -0.5);
	mat44 vox2mm = nifti_mat44_mul(thdr.vox_to_ras, zoomMat);
	tvx_writer w;
	writer_init(&w, hdr, nifti_mat44_mul(nifti_mat44_inverse(sform(hdr)), vox2mm));
	int stride = 3 + thdr.n_scalars;
	float *buf = NULL;
	size_t bufcap = 0;
	for (int i = 0; i < thdr.n_count; i++) {
		int32_t m;
		if (fread(&m, sizeof(m), 1, fp) != 1)
			break;
		if ((size_t)m * stride > bufcap)
			buf = realloc(buf, (bufcap = (size_t)m * stride) * sizeof(float));
		if (fread(buf, sizeof(float), (size_t)m * stride, fp) != (size_t)m * stride)
			break;
		for (int j = 0; j < m; j++)
			add_vertex(&w, buf + j * stride);
		fseek(fp, thdr.n_properties * sizeof(float), SEEK_CUR);
		end_streamline(&w);
	}
	free(buf);
	fclose(fp);
	write_tvx(fnm, &w, delta);
	return EXIT_SUCCESS;
}

static int load_tck(const char *fnm, nifti_1_header *hdr, bool delta) {
	FILE *fp = fopen(fnm, "rb");
	if (fp == NULL)
		return EXIT_FAILURE;
	char line[1024];
	do { // text header ends with "END\n"; assumes Float32LE data follows immediately
		if (!fgets(line, sizeof(line), fp)) {
			fclose(fp);
			return EXIT_FAILURE;
		}
	} while (strcmp(line, "END\n") != 0);
	tvx_writer w;
	writer_init(&w, hdr, nifti_mat44_inverse(sform(hdr)));
	float xyz[3];
	while (fread(xyz, sizeof(xyz), 1, fp) == 1) {
		if (isfinite(xyz[0])) {
			add_vertex(&w, xyz);
			continue;
		}
		end_streamline(&w); // NaN separates streamlines, Inf terminates the file
		if (!isnan(xyz[0]))
			break;
	}
	end_streamline(&w);
	fclose(fp);
	write_tvx(fnm, &w, delta);
	return EXIT_SUCCESS;
}

// ---- query
static bool read_tvx(const char *fnm, tvx_t *t) {
	// gzread transparently reads both raw and gzipped files
	gzFile fgz = gzopen(fnm, "rb");
	if (!fgz || gzread(fgz, &t->h, sizeof(t->h)) != sizeof(t->h) || (t->h.signature != kSigRaw && t->h.signature != kSigDelta)) {
		printf("Not a TVX file: %s\n", fnm);
		return false;
	}
	gzbuffer(fgz, 1 << 20);
	size_t nbytes = t->h.signature == kSigRaw ? (size_t)t->h.nvoxel * 4 : (size_t)t->h.nvoxel * 11;
	t->offsets = malloc(sizeof(uint32_t) * t->h.noffset);
	t->verts = malloc(sizeof(uint32_t) * t->h.nvoxel);
	gzread(fgz, t->offsets, sizeof(uint32_t) * t->h.noffset);
	if (t->h.signature == kSigDelta) {
		uint8_t *buf = malloc(nbytes);
		gzread(fgz, buf, nbytes);
		delta_decode(t, buf);
		free(buf);
	} else
		gzread(fgz, t->verts, nbytes);
	gzclose(fgz);
	return true;
}

static void free_tvx(tvx_t *t) {
	free(t->offsets);
	free(t->verts);
	t->offsets = t->verts = NULL;
}

static bool tvx_matches(const tvx_t *t, const nifti_1_header *hdr) {
	bool ok = true;
	for (int i = 0; i < 3; i++)
		ok &= t->h.dim[i] == (uint32_t)hdr->dim[i + 1];
	for (int i = 0; i < 4; i++)
		ok &= t->h.srow_x[i] == hdr->srow_x[i] && t->h.srow_y[i] == hdr->srow_y[i] && t->h.srow_z[i] == hdr->srow_z[i];
	if (!ok) {
		printf("NIfTI and TVX do not match (use fslhd for NIfTI):\n");
		printf(" dim123: %u %u %u\n", t->h.dim[0], t->h.dim[1], t->h.dim[2]);
		printf(" sto_xyz1: %g %g %g %g\n", t->h.srow_x[0], t->h.srow_x[1], t->h.srow_x[2], t->h.srow_x[3]);
		printf(" sto_xyz2: %g %g %g %g\n", t->h.srow_y[0], t->h.srow_y[1], t->h.srow_y[2], t->h.srow_y[3]);
		printf(" sto_xyz3: %g %g %g %g\n", t->h.srow_z[0], t->h.srow_z[1], t->h.srow_z[2], t->h.srow_z[3]);
	}
	return ok;
}

// fraction of streamlines touching any non-zero voxel of img
static float query_tvx(const tvx_t *t, const uint8_t *img) {
	uint32_t nstreamline = t->h.noffset - 1;
	uint32_t hits = 0;
	for (uint32_t i = 0; i < nstreamline; i++)
		for (uint32_t j = t->offsets[i]; j < t->offsets[i + 1]; j++)
			if (img[t->verts[j]]) {
				hits++;
				break;
			}
	return (float)hits / (float)nstreamline;
}

static void show_help(char *fname) {
	printf("nii2tvx %s %s %s\n", kdate, kCCsuf, kCPUsuf);
	printf("Computes overlap of lesion (NIfTI) and tracts (tvx).\n");
	printf("Usage to create TVX file(s)\n");
	printf(" %s template.nii tracks1.tck tracks2.tck\n", fname);
	printf(" %s -d template.nii tracks1.trk   (-d: delta-encoded TVX, ~4x smaller)\n", fname);
	printf("Usage to compute lesion overlap(s) with TVX file(s)\n");
	printf(" %s lesion.nii tracks1.tvx tracks2.tvx\n", fname);
	printf(" %s lesion1.nii lesion2.nii tracks1.tvx tracks2.tvx\n", fname);
	printf(" %s ./imgs/w*lesion.nii.gz ./tvx/*.tvx > results.tsv\n", fname);
	printf(" -m: low memory, re-read each TVX per lesion instead of keeping all in RAM\n");
	exit(EXIT_FAILURE);
}

static bool is_nifti(const char *fnm) {
	return is_ext(fnm, ".nii") || is_ext(fnm, ".nii.gz");
}

int main(int argc, char **argv) {
	bool delta = false, lowmem = false;
	int first = 1;
	for (; first < argc && argv[first][0] == '-'; first++) {
		if (strcmp(argv[first], "-d") == 0) delta = true;
		else if (strcmp(argv[first], "-m") == 0) lowmem = true;
		else show_help(argv[0]);
	}
	int nnifti = 0, ntrack = 0;
	for (int i = first; i < argc; i++) {
		if (access(argv[i], F_OK) != 0) {
			printf("Unable to find file named '%s'\n", argv[i]);
			exit(EXIT_FAILURE);
		}
		if (is_nifti(argv[i])) nnifti++;
		else ntrack++;
	}
	if (nnifti == 0 || ntrack == 0) {
		printf("Arguments must include at least one NIfTI image and at least one tractography file (TCK, TRK, TVX)\n");
		show_help(argv[0]);
	}
	tvx_t *tvx = calloc(argc, sizeof(tvx_t)); // indexed by argv position, loaded on first use
	float *fracs = malloc(sizeof(float) * argc);
	int *idxs = malloc(sizeof(int) * argc);
	bool header_written = false;
	for (int i = first; i < argc; i++) {
		if (!is_nifti(argv[i]))
			continue;
		nifti_1_header hdr;
		uint8_t *img = load_nii_mask(argv[i], &hdr);
		if (img == NULL) {
			printf("Unable to load NIfTI image %s\n", argv[i]);
			exit(EXIT_FAILURE);
		}
		int nfrac = 0;
		for (int j = first; j < argc; j++) {
			if (is_nifti(argv[j]))
				continue;
			if (is_ext(argv[j], ".tck"))
				load_tck(argv[j], &hdr, delta);
			else if (is_ext(argv[j], ".trk"))
				load_trk(argv[j], &hdr, delta);
			else if (is_ext(argv[j], ".tvx")) {
				if (!tvx[j].verts && !read_tvx(argv[j], &tvx[j]))
					exit(EXIT_FAILURE);
				if (!tvx_matches(&tvx[j], &hdr))
					exit(EXIT_FAILURE);
				fracs[nfrac] = query_tvx(&tvx[j], img);
				idxs[nfrac++] = j;
				if (lowmem)
					free_tvx(&tvx[j]);
			} else {
				printf("Extension unknown %s\n", argv[j]);
				exit(EXIT_FAILURE);
			}
		}
		free(img);
		if (nfrac == 0)
			continue;
		char basenm[768];
		if (!header_written) {
			header_written = true;
			printf("id");
			for (int k = 0; k < nfrac; k++) {
				strcpy(basenm, argv[idxs[k]]);
				strip_ext2(basenm);
				printf("\t%s", basenamex(basenm));
			}
			printf("\n");
		}
		strcpy(basenm, argv[i]);
		strip_ext2(basenm);
		printf("%s", basenamex(basenm));
		for (int k = 0; k < nfrac; k++)
			printf("\t%g", fracs[k]);
		printf("\n");
	}
	for (int j = first; j < argc; j++)
		free_tvx(&tvx[j]);
	free(tvx);
	free(fracs);
	free(idxs);
	exit(EXIT_SUCCESS);
}
