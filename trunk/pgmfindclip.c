/* 
   pgmfindclip.c 

   written by Christian Vogelgsang <chris@lallafa.de>
   under the GNU Public License V2

   modified 2006/09/06 by Nicolas Botti <rududu@laposte.net>
   $Id: pgmfindclip.c,v 1.14 2006/09/06 12:28:05 cnvogelg Exp $
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>

/* safety borders */
static int xsafety = 0;
static int ysafety = 0;

/* verbose mode */
static int verbose = 0;

/* modulos */
static int xfmodulo = 1;
static int yfmodulo = 1;
static int xbmodulo = 1;
static int ybmodulo = 1;

/* expand */
static int expand = 0;

/* luminance only */
static int lumionly = 0;

/* output mplayer crop format */
static int mpformat = 0;

/* search offset */
static int xoffset = 0;
static int yoffset = 0;

/* gnuplot output */
static int gnuplot = 0;
#define GNUPLOT_DATA "temp.dat"
#define GNUPLOT_BIN "gnuplot"

/* Threshold */
#define VAR_THRES 200

/* write pgm */
static int writePGM = 0;

/* clip structure */

typedef struct {
	int t, b, l, r;
} clip;

/* export data set as eps file via gnuplot */
void saveGnuplot(char *name,
		 int *data,int len,
		 int clip0,int clip1,
		 int offset,int thres)
{
  FILE *fh;
  int i;

  /* write data values to file */
  fh = fopen(GNUPLOT_DATA,"w");
  if(fh==0) {
    fprintf(stderr,"Error opening file '%s'\n",GNUPLOT_DATA);
    return;
  }
  for(i=0;i<len;i++)
    fprintf(fh,"%d %d\n",i,data[i]);
  fprintf(fh,"\n\n%d %d\n%d %d\n",clip0,data[clip0],clip1,data[clip1]);
  fclose(fh);

  /* open a pipe to gnuplot command -> generate eps */
  fh = popen(GNUPLOT_BIN,"w");
  fprintf(fh,
	  "set terminal postscript eps\n"

	  "set output \"%s\"\n"
	  "set title \"pgmfindclip plot '%s'\"\n"
	  "plot [%d:%d] "
	  "\"%s\" index 0 title \"gradsum\" with lines lt 1,  "
	  "%d title \"threshold\" with lines lt 2,  "
	  "\"%s\" using 1:(($2>=%d)&&($1>=%d)&&($1<=%d)?$2:1/0) "
	  "title \"candidates\" with impulses lt 1,  "
	  "\"%s\" index 1 title \"clip\" with points pt 6\n",

	  name,
	  name,
	  -10,len+9,
	  GNUPLOT_DATA,
	  thres,
	  GNUPLOT_DATA,thres,offset,len-1-offset,
	  GNUPLOT_DATA);
  fclose(fh);
  unlink(GNUPLOT_DATA);
}

/* ----- calc gradients ----- */

/* apply this filter : [1 -1] */
void calcFilter(int *memory,int w)
{
	int x;

	/* for each pixel in the row */
	for(x=0;x<(w-1);x++) {
		memory[x] = abs(memory[x] - memory[x+1]);
	}
}

/* calc value sum of a row */
void calcStatRow(unsigned char *memory, int w, int y, int * rsum, int * rvar)
{
	int sum = 0;
	long long sum2 = 0;
	int x;
	
	/* current row */
	unsigned char *ptr = memory + y * w;
	/* for each pixel in the row */
	for(x=0;x<w;x++) {
		sum += *ptr;
		sum2 += *ptr * *ptr;
		ptr++;
	}
	*rsum = sum * 100 / w;
	*rvar = sqrtf((float)(sum2 * w - sum * sum)) * 100 / w;
}

/* calc value sum of a column */
void calcStatColumn(unsigned char *memory, int w, int h, int x, int * rsum, int * rvar)
{
	int sum = 0;
	long long sum2 = 0;
	int y;
	
	/* current column */
	unsigned char *ptr = memory + x;
	/* for each pixel in the column */
	for(y=0;y<h;y++) {
		sum += *ptr;
		sum2 += *ptr * *ptr;
		ptr+=w;
	}
	*rsum = sum * 100 / h;
	*rvar = sqrtf((float)(sum2 * h - sum * sum)) * 100 / h;
}

/* Find the clip value */
/* memory : start position */
/* w : lenth to search */
/* dir : direction for finding (1 or -1) */
/* retThres : return Thres used to find the candidate */
int findClipCandidate(int *memory, int * var, int w, int dir, int * retThres)
{
	int i, j = 0, istart, jstart, state = 0, zeros = 0;
	int sum = 0, count = 0, thres;
	long long sum2 = 0;
	for (i = 0; i < w; i++, j += dir)
		if (var[j] > VAR_THRES)
			break;
	if (i==w)
		return -1;
	i -= 4;
	j -= 4 * dir;
	if (i < 0){
		 i = 0;
		 j = 0;
	}
	if (i+32 < w)
		w = i+32;
	istart = i;
	jstart = j;
	for (; i < w; i++, j += dir){
		sum += memory[j];
		sum2 += memory[j] * memory[j];
	}
	sum2 = sqrtf((float)(sum2 * (w-istart) - sum * sum) / ((w-istart) * (w-istart)));
	sum = sum / (w-istart);
	thres = sum + 2 * sum2;

	sum2 = sum = 0;
	for (i = istart, j = jstart; i < w; i++, j += dir){
		if (memory[j] < thres){
			sum += memory[j];
			sum2 += memory[j] * memory[j];
			count++;
		}
	}

	if (count == 0)
		return -1;
	
	sum2 = sqrtf((float)(sum2 * count - sum * sum) / (count * count));
	sum = sum / count;
	thres = 2 * sum + 6 * sum2;
	*retThres = thres;
	
	for (i = istart, j = jstart; i < w; i++, j += dir){
		if (var[j] <= VAR_THRES)
			zeros++;
		if (memory[j] > thres){
			thres = sum + 3 * sum2;
			state |= 1;
		}else{
			if (state & 1 /*&& memory[j+dir] <= thres*/)
				break;
		}
	}
	if (!(state & 1)){
		if (istart == 0)
			return 0;
		return -1;
	}else if (istart == 0 && i > 8 && zeros < 3)
		return 0;
	return i;
}

int findClip(int *memory, int * var, int w, int dir)
{
	int Candidates[8];
	int Thres[8];
	int NbCand = 1, i, j = 0, state=0;
	int MaxThres = 0;

	Candidates[0] = findClipCandidate(memory, var, w, dir, Thres);

	for (i = 0; i < w; i++, j+=dir){
		if (var[j] > VAR_THRES){
			if (var[j] > VAR_THRES * 4)
				break;
			state=1;
		} else if (state > 0){
			state++;
			if (state > 3){
				state = 0;
				Candidates[NbCand] = findClipCandidate(memory + j, var + j, w - i, dir, Thres + NbCand);
				if (Candidates[NbCand] >= 0)
					Candidates[NbCand] += i;
				NbCand++;
				if (NbCand == 8)
					break;
			}
		}
	}

	j = -1;
	for (i = 0; i < NbCand; i++){
		if (Candidates[i] >= 0 && Thres[i] > MaxThres){
			MaxThres = Thres[i];
			j = Candidates[i];
		}
	}
	return j;
}


/* ----- the main clipping algorithm ----- */
int findClipBorders(char *name,
		    unsigned char *memory,int width,int height,
		    int *l,int *r,int *t,int *b)
{
	int *yd, *xd, *yv, *xv;
	int x,y;
	int w,h;
	int fullx = 0;
	int fully = 0;
// 	FILE *fh;
	
	/* reset values */
	*l = *r = width;
	*t = *b = height;
	
	/* number of differences */
	w = width;
	h = height;
	
	/* ydiffs and xdiffs */
	yd = (int *)malloc(sizeof(int)*h);
	xd = (int *)malloc(sizeof(int)*w);
	
	/* variances */
	yv = (int *)malloc(sizeof(int)*h);
	xv = (int *)malloc(sizeof(int)*w);
	
	
	/* calc all sums */
	for(x=0;x < width;x++)
		calcStatColumn(memory,width,height,x, xd + x, xv + x);
	
	for(y=0;y < height;y++)
		calcStatRow(memory,width,y, yd + y, yv + y);
	
	/* calc grad */
	calcFilter(xd, width);
	calcFilter(yd, height);

// 	fh = fopen("h.dat","w");
// 	if(fh==0) {
// 		fprintf(stderr,"Error opening file \n");
// 		return 0;
// 	}
// 	for(x=0;x<width;x++)
// 		fprintf(fh,"%d %d %d\n",x,xd[x],xv[x]);
// 	fclose(fh);
// 
// 	fh = fopen("v.dat","w");
// 	if(fh==0) {
// 		fprintf(stderr,"Error opening file \n");
// 		return 0;
// 	}
// 	for(y=0;y<width;y++)
// 		fprintf(fh,"%d %d %d\n",y,yd[y],yv[y]);
// 	fclose(fh);
	
	/* verbose mode */
	if(verbose) {
		fprintf(stderr,"xsumgrad: ");
		for(x=0;x<w-1;x++)
		fprintf(stderr,"%d ",xd[x]);
		fprintf(stderr,"\nysumgrad: ");
		for(y=0;y<h-1;y++)
		fprintf(stderr,"%d ",yd[y]);
		fprintf(stderr,"\n");
	}
	
	/* ----- vertical operation ----- */
	if (yoffset == 0)
		yoffset = height - 1;
	if (yoffset > height - 1)
		yoffset = height - 1;
	*t = findClip(yd, yv, yoffset, 1);
	*b = findClip(yd + height - 2, yv + height - 1, yoffset, -1);
	
	if (*t < 0 || *b < 0 || (height - *b - *t) < (height >> 2))
		fully = 1;
	
	/* ----- horizontal operation ----- */
	if (xoffset == 0)
		xoffset = width - 1;
	if (xoffset > width - 1)
		xoffset = width - 1;
	*l = findClip(xd, xv, xoffset, 1);
	*r = findClip(xd + width - 2, xv + width - 1, xoffset, -1);
	
	if (*l < 0 || *r < 0 || (width - *r - *l) < (width >> 2))
		fullx = 1;

  /* gnuplot output */
  if(gnuplot) {
    char *dot,*file;
    int len;

    /* construct output name */
    dot = strrchr(name,'.');
    if(dot!=NULL) *dot = '\0';
    len = strlen(name);
    file = (char *)malloc(len+7);
    if(file!=NULL) {
      strcpy(file,name);
      strcpy(file+len,"-x.eps");
      saveGnuplot(file,xd,w,*l,w-1-*r,xoffset,0);
      file[len+1] = 'y';
      saveGnuplot(file,yd,h,*t,h-1-*b,yoffset,0);
      free(file);
    }
  }

  free(xd);
  free(yd);

  return(!(fullx || fully));
}

/* read the next header line of pgm file */
unsigned char *readHeaderLine(FILE *fh)
{
  static unsigned char buf[1024];

  do {
    /* read error */
	  if(fgets((char *)buf,1023,fh)==0) {
      fprintf(stderr,"Read Error!\n");
      return 0;
    }
  } 
  /* skip comment or empty lines */
  while((buf[0]=='#')||(buf[0]=='\n'));

  return buf;
}

/* read a gray level pgm image */
unsigned char *readPGMFile(const char *file,int *width,int *height)
{
	unsigned char *memory;
	unsigned char *line;
	FILE *fh;
	int size;
	int dummy;
	
	/* open pgm image */
	fh = fopen(file,"r");
	if(fh==0) {
		fprintf(stderr,"Can't open '%s'\n",file);
		return 0;
	}
	/* read pgm header */
	line = readHeaderLine(fh);
	if(strcmp("P5\n",(char *)line)!=0) {
		fprintf(stderr,"No P5 pgm!\n");
		return 0;
	}
	/* read image dimension */
	line = readHeaderLine(fh);
	if(sscanf((char *)line,"%d %d %d",width,height,&dummy)!=3) {
		if(sscanf((char *)line,"%d %d",width,height)!=2) {
			fprintf(stderr,"No dimension in pgm found!\n");
			fclose(fh);
			return 0;
		}
		/* skip component max value */
		line = readHeaderLine(fh);
	}
	
	/* load only the lumi part of a mplayer pgm file */
	if(lumionly)
		*height = *height * 2 / 3;
	
	/* get memory */
	size = *width * *height;
	memory = (unsigned char *)malloc(size);
	if(memory==0) {
		fprintf(stderr,"Out of memory!\n");
		fclose(fh);
		return 0;
	}
	/* read data */
	if(fread(memory,size,1,fh)!=1) {
		fprintf(stderr,"Read error while fetching raw data!\n");
		free(memory);
		fclose(fh);
		return 0;
	}
	
	/* close file */
	fclose(fh);
	
	return memory;
}

/* fmodulo alignment - align to smaller block */
int modalign(int value,int fmodulo,int offset)
{
  if (fmodulo != 0)
    {
      int block = (value+offset)/fmodulo;
      return block * fmodulo;
    } else
      return value+offset;
}

/* align the frame and border */
void alignFrame(int *l,int *r,int w,int fmod,int bmod)
{
  int nw,offset,aw,border,blocks,left,right;

  /* new exact width */
  nw = w - (*l + *r);
  /* align width */
  offset = expand ? fmod - 1 : 0;
  aw = modalign(nw,fmod,offset);
  if(aw>w) aw -= fmod;

  /* calc number of border pixels */
  border = w - aw;
  if(border==0) {
    *l = *r = 0;
    return;
  }

  /* check if border alignment is possible */
  if((border % bmod)!=0) {
    fprintf(stderr,"Warning: frame AND border alignment impossible!\n"
	    " -> Ignoring border alignment!\n");
    bmod = 1;
  }

  /* calc border blocks */
  blocks = border / bmod;

  /* distribute blocks with the ratio of the exact borders */
  if (*l + *r > 0) {
    left  = blocks * *l / (*l + *r);
    right = blocks - left;
    left  *= bmod;
    right *= bmod;
  }
  else left = right = 0;

  /* check values */
  if(expand) {
    if((*l<left)&&(*r<right)) {
      fprintf(stderr,"Warning: no real expansion!\n");
    }
  } else {
    if((*l>left)||(*r>right)) {
      fprintf(stderr,"Warning: no real shrink!\n");
    }
  }

  *l = left;
  *r = right;
}

int clipComp(clip * a, clip * b)
{
	if (a->t == b->t && a->b == b->b && a->l == b->l && a->r == b->r)
		return 1;
	return 0;
}

/* filters all the valid borders found */
void findBorders(clip * clips, int count)
{
	int * counts = (int*)malloc(sizeof(int) * count);
	clip * b = (clip*)malloc(sizeof(clip) * count);
	clip sum = {0,0,0,0};
	clip sum2 = {0,0,0,0};
	int i, c;
	int min_count;

	if (count == 1)
		return;

	b[0] = clips[0];
	c = 1;
	counts[c] = 1;

	for (i = 0; i < count; i++){
		int j = 0;
		for (j = 0; j < c; j++)
			if (clipComp(b+j,clips+i))
				break;
		if (j == c){
			b[c] = clips[i];
			counts[c] = 1;
			c++;
		}else{
			counts[j]++;
		}
	}

	min_count = count / 64;
	count = 0;

	for (i = 0; i < c; i++){
		if (counts[i] >= min_count){
			counts[count] = counts[i];
			b[count] = b[i];
			count++;
		}
	}

	c = count;
	count = 0;

	if(verbose)
		for (i = 0; i < c; i++)
			fprintf(stderr,"[%d,%d,%d,%d] : %d\n",b[i].t,b[i].l,b[i].b,b[i].r,counts[i]);

	for (i = 0; i < c; i++){
		sum.t += b[i].t * counts[i];
		sum.b += b[i].b * counts[i];
		sum.l += b[i].l * counts[i];
		sum.r += b[i].r * counts[i];
		sum2.t += b[i].t * b[i].t * counts[i];
		sum2.b += b[i].b * b[i].b * counts[i];
		sum2.l += b[i].l * b[i].l * counts[i];
		sum2.r += b[i].r * b[i].r * counts[i];
		count += counts[i];
	}

	/* standard deviation * 100 */
	sum2.t = sqrtf((float)((float)sum2.t * count - sum.t * sum.t) / (count * count)) * 100;
	sum2.b = sqrtf((float)((float)sum2.b * count - sum.b * sum.b) / (count * count)) * 100;
	sum2.l = sqrtf((float)((float)sum2.l * count - sum.l * sum.l) / (count * count)) * 100;
	sum2.r = sqrtf((float)((float)sum2.r * count - sum.r * sum.r) / (count * count)) * 100;
	/* valid value range */
	sum.t = (sum.t * 100 / count + sum2.t * 2 + 50) / 100;
	sum.b = (sum.b * 100 / count + sum2.b * 2 + 50) / 100;
	sum.l = (sum.l * 100 / count + sum2.l * 2 + 50) / 100;
	sum.r = (sum.r * 100 / count + sum2.r * 2 + 50) / 100;

	sum2.t = sum2.b = sum2.l = sum2.r = 0;

	for (i = 0; i < c; i++){
		if (b[i].t > sum2.t && b[i].t <= sum.t)
			sum2.t = b[i].t;
		if (b[i].b > sum2.b && b[i].b <= sum.b)
			sum2.b = b[i].b;
		if (b[i].l > sum2.l && b[i].l <= sum.l)
			sum2.l = b[i].l;
		if (b[i].r > sum2.r && b[i].r <= sum.r)
			sum2.r = b[i].r;
	}

	clips[0].t =sum2.t;
	clips[0].b =sum2.b;
	clips[0].l =sum2.l;
	clips[0].r =sum2.r;

	free (counts);
	free (b);
}

/* move and align border values */
void modifyBorders(int *l,int *r,int *t,int *b,int width,int height)
{
  if(verbose)
    fprintf(stderr,"found:   %d,%d,%d,%d\n",*t,*l,*b,*r);

  /* add safety borders */
  *t += ysafety;
  *b += ysafety;
  *l += xsafety;
  *r += xsafety;

  if(verbose)
    fprintf(stderr,"safety:  %d,%d,%d,%d\n",*t,*l,*b,*r);

  if(verbose)
    fprintf(stderr,"x-align\n");

  /* do frame/border alignment */
  if(xfmodulo!=1)
    alignFrame(l,r,width,xfmodulo,xbmodulo);
  /* only border alignment */
  else if(xbmodulo!=1) {
    int offset = expand ? 0 : xbmodulo-1;
    *l = modalign(*l,xbmodulo,offset);
    *r = modalign(*r,xbmodulo,offset);
  }

  if(verbose)
    fprintf(stderr,"y-align\n");

  /* do frame/border alignment */
  if(yfmodulo!=1)
    alignFrame(t,b,height,yfmodulo,ybmodulo);
  /* only border alignment */
  else if(ybmodulo!=1) {
    int offset = expand ? 0 : ybmodulo-1;
    *t = modalign(*t,ybmodulo,offset);
    *b = modalign(*b,ybmodulo,offset);
  }

  if(verbose)
    fprintf(stderr,"aligned: %d,%d,%d,%d\n",*t,*l,*b,*r);
}

/* write border markers into image */
void writeMarkers(char *name,unsigned char *memory,int w,int h,
		  int l,int r,int t,int b)
{
  char *dot,*file;
  int len;
  int x,y;
  char *ptr1,*ptr2;

  /* draw vertical markers */
  ptr1 = (char *)memory +  t    * w + l;
  ptr2 = (char *)memory + (t+1) * w - 1 - r;
  for(y=t;y<(h-b);y++) {
    *ptr1 = 255 - *ptr1;
    *ptr2 = 255 - *ptr2;
    ptr1 += w;
    ptr2 += w;
  }

  /* draw horizontal markers */
  ptr1 = (char *)memory +  t      * w + l + 1;
  ptr2 = (char *)memory + (h-1-b) * w + l + 1;
  for(x=l+1;x<(w-r-1);x++) {
    *ptr1 = 255 - *ptr1;
    *ptr2 = 255 - *ptr2;
    ptr1++;
    ptr2++;
  }

  /* construct output name */
  dot = strrchr(name,'.');
  if(dot!=NULL) *dot = '\0';
  len = strlen(name);
  file = (char *)malloc(len+7);
  if(file!=NULL) {
    FILE *fh;

    strcpy(file,name);
    strcpy(file+len,"-m.pgm");

    /* save pgm file */
    fh = fopen(file,"w");
    if(fh!=NULL) {
      fprintf(fh,"P5\n%d %d 255\n",w,h);
      fwrite(memory,w*h,1,fh);
      fclose(fh);
    }
  }
}

/* print usage and exit */
void usage(void)
{
  fprintf(stderr,
"pgmfindclip - find clipping borders\n"
"  written by Christian Vogelgsang <chris@lallafa.de>\n"
"  under the GNU Public License V2\n"
"  $Revision: 1.13 $\n\n"	  
"Usage: pgmfindclip [Options] <image.pgm> ...\n\n"
" -s <safety>[,<ysafety>]   add a safety border     (default: none)\n"
" -b <modulo>[,<ymodulo>]   align clip borders      (default: none)\n"
" -f <modulo>[,<ymodulo>]   align output frame size (default: none)\n"
" -o <offset>[,<yoffset>]   search begin offset     (default: 0 0)\n"
" -e                        expand and do not shrink in frame alignment\n"
" -y                        input pgm files are yuv files from mplayer\n"
" -m                        output formated for mencoder crop\n"
" -v                        enable verbose mode\n"
" -p                        plot results to *x.eps, *y.eps (req. gnuplot)\n"
" -w                        write PGM files with border markers (*-m.pgm)\n"
		 );
  exit(1);
}

/* parse an optional two value integer argument */
void twoValArg(char *str,int *v1,int *v2)
{
  char *comma = strchr(str,',');
  if(comma==NULL) {
    *v1 = *v2 = atoi(str);
  } else {
    *comma = '\0';
    comma++;
    *v1 = atoi(str);
    *v2 = atoi(comma);
  }
}

/* loop over all images */
int main(int argc,char *argv[])
{
  int i;
  int width,height,left,right,top,bottom;
  int count;
  clip * clips;

  /* parse args */
  int c;
  while((c=getopt(argc,argv,"s:f:b:o:veypwm"))!=-1) {
    switch(c) {
    case 's':
      twoValArg(optarg,&xsafety,&ysafety);
      break;
    case 'f':
      twoValArg(optarg,&xfmodulo,&yfmodulo);
      break;
    case 'b':
      twoValArg(optarg,&xbmodulo,&ybmodulo);
      break;
    case 'o':
      twoValArg(optarg,&xoffset,&yoffset);
      break;
    case 'v':
      verbose = 1;
      break;
    case 'e':
      expand = 1;
      break;
    case 'y':
      lumionly = 1;
      break;
	case 'm':
	  mpformat = 1;
	  break;
    case 'p':
      gnuplot = 1;
      break;
    case 'w':
      writePGM = 1;
      break;
    default:
      usage();
    }
  }

  /* no images given */
  if(optind == argc)
    usage();

  /* loop over images */
  count = 0;
  clips = (clip *)malloc(sizeof(clip)*(argc-optind));
  for(i=optind;i<argc;i++) {
    int l,r,t,b;
    unsigned char *memory;

    /* read image */
    memory = readPGMFile(argv[i],&width,&height);
    if(memory==0) {
      fprintf(stderr,"Error reading PGM image '%s'\n",argv[i]);
      exit(1);
    }

    /* ----- Find Crop Border ----- */
    if(findClipBorders(argv[i],memory,width,height,&l,&r,&t,&b)) {
      if(verbose)
		fprintf(stderr,"image (%d x %d) clip: t=%d l=%d b=%d r=%d\n",
			width,height,t,l,b,r);

	  clips[count].t = t;
	  clips[count].b = b;
	  clips[count].l = l;
	  clips[count].r = r;

      /* count valid images */
      count++;
    } else {
      if(verbose)
		fprintf(stderr,"no clip region found!\n");
    }

    /* free memory */
    free(memory);
  }
  if(count==0) {
    fprintf(stderr,"No valid images found!\n");
    exit(1);
  }

  findBorders(clips, count);
  
  left = clips[0].l;
  right = clips[0].r;
  top = clips[0].t;
  bottom = clips[0].b;
  
  free(clips);
  /* move and align the clip borders */
  modifyBorders(&left,&right,&top,&bottom,width,height);

  /* write PGM images with clip borders */
  if(writePGM) {
    for(i=optind;i<argc;i++) {
      unsigned char *memory = readPGMFile(argv[i],&width,&height);
      writeMarkers(argv[i],memory,width,height,left,right,top,bottom);
      free(memory);
    }
  }

  if (mpformat)
	  printf("%d:%d:%d:%d\n",width - left - right, height - top - bottom, left,top);
  else
	/* generate transcode-friendly result */
	printf("%d,%d,%d,%d\n",top,left,bottom,right);


  exit(0);
}
