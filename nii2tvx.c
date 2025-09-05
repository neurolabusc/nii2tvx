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

void change_ext(char* input, char* new_extension, int size) {
	char* output = input; // save pointer to input in case we need to append a dot and add at the end of input
	while(*(++input) != '\0') // move pointer to final position
		;
	while(*(--input) != '.' && --size > 0) // start going backwards until we encounter a dot or we go back to the start
		;
	// if we've encountered a dot, let's replace the extension, otherwise let's append it to the original string
	size == 0 ? strncat(output, new_extension, 4 ) : strncpy(input, new_extension, 4);
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

char *basenamex(char const *path) {
	char *s = strrchr(path, '/');
	if(s==NULL)
		return strdup(path);
	else
		return strdup(s + 1);
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

#define MY_LITTLE_ENDIAN (((union { unsigned x; unsigned char c; }){1}).c)
bool isnanx(float f) { //isnan disabled by gcc -Ofast and -ffinite-math-only
//4byte IEEE: msb[31] = signbit, bits[23-30] exponent, bits[0..22] mantissa
//exponent of all 1s = Infinity, NAN or Indeterminate
	uint32_t i = *(long *)&f;
	#ifdef MY_LITTLE_ENDIAN
		return ((i&0x7f800000)==0x7f800000)&&(i&0x007fffff);
	#else
		return ((i&0x0000807f)==0x0000807f)&&(i&0xffff7f00);
	#endif
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
	uint32_t signature; //must be 175666804 "TVX\n"
	uint32_t dim[3]; //correspond to NIfTI dim[1..3]
	float srow_x[4] ; // 1st row affine transform.
	float srow_y[4] ; // 2nd row affine transform.
	float srow_z[4] ; // 3rd row affine transform.
	uint32_t noffset; //number of offsets (streamlines +1)
	uint32_t nvoxel; //number of voxels stored for all streamlines
};
#pragma pack()
typedef struct tvx_header tvx_header;

/*void reportMat(mat44 m) {
	printf("m = [");
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			printf(" %g", m.m[i][j]);
			if ((j == 3) && (i < 3))
				printf(";");
		}
	}
	printf("]\n");
}*/

int write_tvx(const char *fnm, nifti_1_header * hdr, size_t noffset, size_t nvert, uint32_t* offsets, uint32_t* verts) {
	char *outnm = malloc(strlen(fnm) + 1); 
	strcpy(outnm, fnm);
	change_ext(outnm, ".tvx", sizeof(outnm));
	FILE *fp = fopen(outnm,"wb");
	free(outnm);
	//create tvx header
	tvx_header thdr;
	thdr.signature = 175666804;
	thdr.dim[0] = hdr->dim[1];
	thdr.dim[1] = hdr->dim[2];
	thdr.dim[2] = hdr->dim[3];
	for (int i = 0; i < 4; i++) { //multiply Pcrs * m
		thdr.srow_x[i] = hdr->srow_x[i];
		thdr.srow_y[i] = hdr->srow_y[i];
		thdr.srow_z[i] = hdr->srow_z[i];
	}
	thdr.noffset = noffset;
	thdr.nvoxel = nvert;
	//write tvx header
	fwrite(&thdr, sizeof(thdr), 1, fp);
	//for (int i = 0; i < 4; i++)
	//	printf("%d\n",offsets[i] );
	//write offsets
	fwrite(offsets, sizeof(uint32_t) * noffset, 1, fp);
	fwrite(verts, sizeof(uint32_t) * nvert, 1, fp);
	fclose(fp);
	return EXIT_SUCCESS;
} // write_tvx()

int load_trk(const char *fnm, nifti_1_header * hdr) {
	//https://trackvis.org/docs/?subsect=fileformat
	//https://brain.labsolver.org/hcp_trk_atlas.html
	//port of https://github.com/tee-ar-ex/trx-javascript
	FILE *fp = fopen(fnm,"rb");
	if (fp == NULL)
		return EXIT_FAILURE;
	fseek(fp, 0L, SEEK_END);
	size_t file_len = ftell(fp);
	fseek(fp, 0, SEEK_SET);
	trk_header thdr;
	fread(&thdr, sizeof(trk_header), 1, fp);
	if ((thdr.n_count == 0) || (thdr.hdr_size != sizeof(trk_header)) || (thdr.version != 2)) {
		printf("Unable to read TRK header %d %d\n", thdr.hdr_size, thdr.version);
		if (thdr.n_count == 0)
			printf("unable to read TRK with implicit n_count (hint: convert to TCK with tff_convert_tractogram.py)");
		fclose(fp);
		return EXIT_FAILURE;
	}
	mat44 zoomMat = make_mat44(
		1.0/thdr.voxel_size[0], 0, 0, -0.5,
		0, 1.0/thdr.voxel_size[1], 0, -0.5,
		0, 0, 1.0/thdr.voxel_size[2], -0.5);
	//if (false)
	//zoomMat = make_mat44(1,0,0,0, 0,1,0,0, 0,0,1,0);
	if (thdr.vox_to_ras.m[3][3] == 0.0) {
		printf("TRK vox_to_ras not set");
		fclose(fp);
		return EXIT_FAILURE;
	}
	mat44 vox2mmMat = nifti_mat44_mul(thdr.vox_to_ras, zoomMat );
	//reportMat(vox2mmMat);
	mat44 mat = sform(hdr);
	mat44 inv = nifti_mat44_inverse(mat);
	float xyz[3];
	float xyzMM[3];
	int prev_vxl = -1;
	size_t noffset = 0;
	size_t nvert_in_streamline = 0;
	size_t nvert = 0;
	size_t unique_nvert = 0;
	size_t scalar_skip = thdr.n_scalars * sizeof(float);
	size_t property_skip = thdr.n_properties * sizeof(float);
	//overprovision vertices
	uint32_t* verts = (uint32_t*)malloc(file_len - ftell(fp));
	uint32_t* offsets = (uint32_t*)malloc((thdr.n_count + 1) * sizeof(uint32_t));
	offsets[noffset] = 0; //start of first streamline
	noffset++;
	for (int i = 0; i < thdr.n_count; i++) {
		int32_t m;
		fread(&m, sizeof(int32_t), 1, fp);
		for (int j = 0; j < m; j++) {
			fread((&xyz), sizeof(xyz), 1, fp);
			size_t pos = ftell(fp);
			fseek(fp, pos + scalar_skip, SEEK_SET);
			xyzMM[0] = xyz[0]*vox2mmMat.m[0][0] + xyz[1]*vox2mmMat.m[0][1] + xyz[2]*vox2mmMat.m[0][2] + vox2mmMat.m[0][3];
			xyzMM[1] = xyz[0]*vox2mmMat.m[1][0] + xyz[1]*vox2mmMat.m[1][1] + xyz[2]*vox2mmMat.m[1][2] + vox2mmMat.m[1][3];
			xyzMM[2] = xyz[0]*vox2mmMat.m[2][0] + xyz[1]*vox2mmMat.m[2][1] + xyz[2]*vox2mmMat.m[2][2] + vox2mmMat.m[2][3];
			//if (j < 3)
			//	printf(">>%d %d %g %g %g -> %g %g %g\n", m, j, xyz[0], xyz[1], xyz[2], xyzMM[0], xyzMM[1], xyzMM[2]);
			vec3i v3 = nifti_xyz_mat44_mul(xyzMM[0], xyzMM[1], xyzMM[2], inv);
			nvert++;
		if ((v3.v[0] < 0) || (v3.v[1] < 0) || (v3.v[2] < 0))
			continue;
		if ((v3.v[0] >= hdr->dim[1]) || (v3.v[1] >= hdr->dim[2]) || (v3.v[2] >= hdr->dim[3]))
			continue;
		int vxl = v3.v[0] + (v3.v[1] * hdr->dim[1]) + (v3.v[2] * hdr->dim[1] * hdr->dim[2]);
		if (prev_vxl == vxl)
			continue;
		prev_vxl = vxl;
		verts[unique_nvert] = vxl;
		unique_nvert ++;
		nvert_in_streamline ++;
		}
		size_t pos = ftell(fp);
		fseek(fp, pos + property_skip, SEEK_SET);
		if (nvert_in_streamline > 0) {
			//start of next streamline
			offsets[noffset] = unique_nvert;
			noffset++;
		}
		nvert_in_streamline = 0;
	}
	fclose(fp);
	//write offsets
	write_tvx(fnm, hdr, noffset, unique_nvert, offsets, verts);
	//write_tvx(const char *fnm, nifti_1_header * hdr, size_t noffset, size_t nvert, uint32_t* offsets, uint32_t* verts) {
	printf("TRK2TVX\t%zu\tstreamlines\tvertices\t%zu\t(\t%zu\tunique\t)\n", noffset - 1, nvert, unique_nvert);
	free(offsets);
	free(verts);
	return EXIT_SUCCESS;
} // load_trk()

int load_tck(const char *fnm, nifti_1_header * hdr) {
	FILE *fp = fopen(fnm,"rb");
	if (fp == NULL)
		return EXIT_FAILURE;
	//skip variable length header
	size_t header_len = 0;
	char search[] = "END\n";
	size_t search_len = strlen(search);
	char buffer[1024];
	size_t offset = 0;
	size_t bytesRead;
	mat44 mat = sform(hdr);
	mat44 inv = nifti_mat44_inverse(mat);
	while ((header_len == 0) && ((bytesRead = fread(buffer, 1, sizeof(buffer), fp)) > 0)) {
		for (size_t i = 0; i < bytesRead; ++i) {
			if (buffer[i] == search[0]) {
				// Check if the rest of the pattern matches
				if (memcmp(buffer + i, search, search_len) == 0) {
					header_len = offset + i + search_len;
					break;
				}
			}
		}
		offset += bytesRead - search_len;
	}
	if (header_len == 0) {
		fclose(fp);
		return EXIT_FAILURE;
	}
	printf("Header length %zu\n", header_len);
	fseek(fp, 0L, SEEK_END);
	size_t file_len = ftell(fp);
	fseek(fp, header_len, SEEK_SET);
	size_t file_pos = header_len;
	float xyz[3];
	size_t noffset = 0;
	size_t nvert_in_streamline = 0;
	size_t nvert = 0;
	size_t unique_nvert = 0;
	int prev_vxl = -1;
	//overprovision arrays
	uint32_t* offsets = (uint32_t*)malloc(file_len - header_len);
	uint32_t* verts = (uint32_t*)malloc(file_len - header_len);
	offsets[noffset] = 0; //start of first streamline
	noffset++;
	while ((file_pos + sizeof(xyz)) <= file_len) {
		fread((&xyz), sizeof(xyz), 1, fp);
		//printf(">> %g %g %g\n", xyz[0],xyz[1],xyz[2]);
		file_pos += sizeof(xyz);
		//if (isnanx(xyz[0])) {
		if (!isfinite(xyz[0])) {
			if (nvert_in_streamline > 0) {
				//start of next streamline
				offsets[noffset] = unique_nvert;
				noffset++;
			}
			nvert_in_streamline = 0;
			prev_vxl = -1;
			if (!isnanx(xyz[0]))
				break;
		} else {
			nvert++;
			vec3i v3 = nifti_xyz_mat44_mul(xyz[0], xyz[1], xyz[2], inv);
			if ((v3.v[0] < 0) || (v3.v[1] < 0) || (v3.v[2] < 0))
				continue;
			if ((v3.v[0] >= hdr->dim[1]) || (v3.v[1] >= hdr->dim[2]) || (v3.v[2] >= hdr->dim[3]))
				continue;
			int vxl = v3.v[0] + (v3.v[1] * hdr->dim[1]) + (v3.v[2] * hdr->dim[1] * hdr->dim[2]);
			if (prev_vxl == vxl)
				continue;
			prev_vxl = vxl;
			verts[unique_nvert] = vxl;
			unique_nvert ++;
			nvert_in_streamline ++;
			if (nvert < 4)
				printf("%d: %g %g %g -> %d %d %d\n", vxl, xyz[0],xyz[1],xyz[2], v3.v[0], v3.v[1], v3.v[2]);
		}
	}
	printf("TCK %zu streamlines, %zu vertices (%zu unique)\n", noffset-1, nvert, unique_nvert);
	fclose(fp);
	write_tvx(fnm, hdr, noffset, unique_nvert, offsets, verts);
	free(offsets);
	free(verts);
	return EXIT_SUCCESS;
} // load_tck()

float load_tvx(const char *fnm, nifti_1_header * hdr, uint8_t * img) {
	FILE *fp = fopen(fnm,"rb");
	tvx_header thdr;
	fread(&thdr, sizeof(tvx_header), 1, fp);
	const int kGzSig = 0x8B1F; //GZ files start 0x1F8B swap for little-endian
	bool isGz = (thdr.signature & kGzSig) == kGzSig;
	uint32_t* offsets = NULL;
	uint32_t* verts = NULL;
	if (!isGz) {
		if (thdr.signature != 175666804) {
			printf("raw TVX signature is wrong");
			fclose(fp);
			return NAN;
		}
		offsets = (uint32_t*)malloc(sizeof(uint32_t) * thdr.noffset);
		verts = (uint32_t*)malloc(sizeof(uint32_t) * thdr.nvoxel);
		fread(offsets, sizeof(uint32_t) * thdr.noffset, 1, fp);
		fread(verts, sizeof(uint32_t) * thdr.nvoxel, 1, fp);
	}
	fclose(fp);
	if (isGz) {
		gzFile fgz = gzopen(fnm, "r");
		gzread(fgz, &thdr, sizeof(tvx_header));
		if (thdr.signature != 175666804) {
			printf("gzipped TVX signature is wrong");
			gzclose(fgz);
			return NAN;
		}
		offsets = (uint32_t*)malloc(sizeof(uint32_t) * thdr.noffset);
		verts = (uint32_t*)malloc(sizeof(uint32_t) * thdr.nvoxel);
		gzread(fgz, offsets, sizeof(uint32_t) * thdr.noffset);
		gzread(fgz, verts, sizeof(uint32_t) * thdr.nvoxel);
		gzclose(fgz);
	}
	bool is_match = true;
	for (int i = 0; i < 3; i++)
		if (thdr.dim[i] != hdr->dim[i+1])
			is_match = false;
	for (int i = 0; i < 4; i++) {
		if (thdr.srow_x[i] != hdr->srow_x[i])
			is_match = false;
		if (thdr.srow_y[i] != hdr->srow_y[i])
			is_match = false;
		if (thdr.srow_z[i] != hdr->srow_z[i])
			is_match = false;
	}
	if (!is_match) {
		printf("NIfTI and TVX do not match (use fslhd for NIfTI):\n");
		printf(" dim123: %d %d %d \n", thdr.dim[0], thdr.dim[1], thdr.dim[2]);
		printf(" sto_xyz1: %g %g %g %g\n", thdr.srow_x[0], thdr.srow_x[1], thdr.srow_x[2], thdr.srow_x[3]);
		printf(" sto_xyz2: %g %g %g %g\n", thdr.srow_y[0], thdr.srow_y[1], thdr.srow_y[2], thdr.srow_y[3]);
		printf(" sto_xyz3: %g %g %g %g\n", thdr.srow_z[0], thdr.srow_z[1], thdr.srow_z[2], thdr.srow_z[3]);
		return NAN;
	}
	uint32_t nstreamline = thdr.noffset - 1; //-1 for fencepost
	uint32_t nstreamlineInMask = 0;
	for (int i = 0; i < nstreamline; i++) {
		uint32_t startOffset = offsets[i];
		uint32_t endOffset = offsets[i+1] -1;
		if ((int64_t)startOffset >= (int64_t)endOffset)
			continue;
		for (int j = startOffset; j <= endOffset; j++) {
			if (img[verts[j]] != 0) {
				nstreamlineInMask ++;
				break;
			}
		}
	}
	//printf("%d/%d\n", nstreamlineInMask, nstreamline);
	free(offsets);
	free(verts);
	return (float)nstreamlineInMask/ (float) nstreamline;
}

void show_help(char *fname){
		printf("nii2tvk %s %s %s\n", kdate, kCCsuf, kCPUsuf);
		printf("Computes overlap of lesion (NIfTI) and tracts (tvx).\n");
		printf("Supply NIfTI image(s) and tracks to process\n");
		printf("Usage to create TVX file(s)\n");
		printf(" %s template.nii tracks1.tck tracks2.tck\n",fname);
		printf(" %s template.nii tracks1.trk\n",fname);
		printf("Usage to compute lesion overlap(s) with TVX file(s)\n");
		printf(" %s lesion.nii tracks1.tvx tracks2.tvx\n",fname);
		printf(" %s lesion1.nii lesion2.nii tracks1.tvx tracks2.tvx\n",fname);
		printf(" %s ./imgs/w*lesion.nii.gz ./tvx/*.tvx > results.tsv\n",fname);
		exit(EXIT_FAILURE);
}

bool is_nifti(char *fnm) {
	return is_ext(fnm,".nii") || is_ext(fnm,".nii.gz");
}

int main(int argc,char **argv) {
	// Check the command line, minimal is name of input and output files
	if ((argc < 2) || (argv[1][0] == '-'))
		show_help(argv[0]);
	int nnifti = 0;
	int nnotnifti = 0;
	for (int i = 1; i < argc; i++) {
		if (!is_nifti(argv[i]))
			continue;
		nnifti++;
		nifti_1_header hdr;
		uint8_t * img = load_nii_mask(argv[i], &hdr);
		if (img == NULL) {
			printf("Unable to load NIfTI image %s\n", argv[i]);
			exit(EXIT_FAILURE);
		}
		float* fracs = (float*)malloc(sizeof(float) * argc);
		int32_t* idxs = (int32_t*)malloc(sizeof(int32_t) * argc);
		int nfrac = 0;
		//apply NIfTI to all non-NIfTI tracks (trx, tck, tvx)
		for (int j = 1; j < argc; j++) {
			if (is_nifti(argv[j]))
				continue;
			if (access(argv[j], F_OK) != 0) {
				printf("Unable to find file named '%s'\n", argv[j]);
				exit(EXIT_FAILURE);
			}
			nnotnifti++;
			if (is_ext(argv[j],".tck"))
				load_tck(argv[j], &hdr);
			else if (is_ext(argv[j],".trk"))
				load_trk(argv[j], &hdr);
			else if (is_ext(argv[j],".tvx")) {
				float frac = load_tvx(argv[j], &hdr, img);
				if (isnan(frac))
					exit(EXIT_FAILURE);
				fracs[nfrac] = frac;
				idxs[nfrac] = j;
				nfrac++;
			} else {
				printf("Extension unknown %s\n", argv[j]);
				exit(EXIT_FAILURE);
			}
		}	
		if (nfrac > 0) {
			char basenm[768];
			if (nnifti == 1) {
				//write header
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
		free(fracs);
		free(idxs);
		free(img);
	}
	if ((nnifti == 0) || (nnotnifti == 0)) {
		printf("Arguments must include at least one NIfTI image and at least one Tractography file (TCK, TRK)");
		show_help(argv[0]);
	}
	exit(EXIT_SUCCESS);
}
