/* 
   pgmfindclip.c 

   written by Christian Vogelgsang <chris@lallafa.de>
   under the GNU Public License V2

   $Id: pgmfindclip.c,v 1.13 2003/05/18 12:28:05 cnvogelg Exp $
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

/* threshold values */
static int xthres = 100;
static int ythres = 100;

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

/* search offset */
static int xoffset = 0;
static int yoffset = 0;

/* gnuplot output */
static int gnuplot = 0;
#define GNUPLOT_DATA "temp.dat"
#define GNUPLOT_BIN "gnuplot"

/* write pgm */
static int writePGM = 0;

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

/* calc difference between row y and y+1 and sum absolute diffs */
int calcSumRowGrad(unsigned char *memory,int w,int y)
{
  int sum = 0;
  int x;

  /* current row */
  unsigned char *ptr = memory + y * w;
  /* next row */
  unsigned char *next = ptr + w;

  /* for each pixel in the row */
  for(x=0;x<w;x++) {
    int diff = abs(*ptr - *next);
    sum += diff;

    /* next column */
    ptr++;
    next++;
  }
  return sum;
}

/* calc difference between cell x and x+1 and sum absolute diffs */
int calcSumColumnGrad(unsigned char *memory,int w,int h,int x)
{
  int sum = 0;
  int y;

  /* current column */
  unsigned char *ptr = memory + x;
  /* next column */
  unsigned char *next = ptr + 1;

  /* for each pixel in the column */
  for(y=0;y<h;y++) {
    int diff = abs(*ptr - *next);
    sum += diff;

    /* next row */
    ptr+=w;
    next+=w;
  }
  return sum;
}

/* calc value sum of a row - grad for border case */
int calcSumRow(unsigned char *memory,int w,int y)
{
  int sum = 0;
  int x;

  /* current row */
  unsigned char *ptr = memory + y * w;
  /* for each pixel in the row */
  for(x=0;x<w;x++) {
    sum += *ptr;
    ptr++;
  }
  return sum;
}

/* calc value sum of a column - grad for border case */
int calcSumColumn(unsigned char *memory,int w,int h,int x)
{
  int sum = 0;
  int y;

  /* current column */
  unsigned char *ptr = memory + x;
  /* for each pixel in the column */
  for(y=0;y<h;y++) {
    sum += *ptr;
    ptr+=w;
  }
  return sum;
}

/* normalize sum and scale */
int norm(int sum,int len)
{
  return(sum * 100 / len);
}

/* ----- the main clipping algorithm ----- */
int findClipBorders(char *name,
		    unsigned char *memory,int width,int height,
		    int *l,int *r,int *t,int *b)
{
  int *yd;
  int *xd;
  int x,y;
  int w,h;
  int fullx;
  int fully;

  /* reset values */
  *l = *r = width;
  *t = *b = height;

  /* number of differences "between" the columns/rows + 2 borders */
  w = width+1;
  h = height+1;

  /* ydiffs and xdiffs */
  yd = (int *)malloc(sizeof(int)*h);
  xd = (int *)malloc(sizeof(int)*w);

  /* calc all differences */
  xd[0] = norm(calcSumColumn(memory,width,height,0),height);
  for(x=0;x<(width-1);x++) 
    xd[x+1] = norm(calcSumColumnGrad(memory,width,height,x),height);
  xd[width] = norm(calcSumColumn(memory,width,height,width-1),height);

  yd[0] = norm(calcSumRow(memory,width,0),width);
  for(y=0;y<(height-1);y++)
    yd[y+1] = norm(calcSumRowGrad(memory,width,y),width);
  yd[height] = norm(calcSumRow(memory,width,height-1),width);

  /* verbose mode */
  if(verbose) {
    fprintf(stderr,"xsumgrad: ");
    for(x=0;x<w;x++)
      fprintf(stderr,"%d ",xd[x]);
    fprintf(stderr,"- xthreshold=%d\n",xthres);
    fprintf(stderr,"ysumgrad: ");
    for(y=0;y<h;y++)
      fprintf(stderr,"%d ",yd[y]);
    fprintf(stderr,"- ythreshold=%d\n",ythres);
  }

  /* ----- vertical operation ----- */
  /* find top clip */
  for(y=yoffset;y<h-yoffset;y++) {
    if(yd[y]>=ythres) {
      *t = y;
      break;
    }
  }
  /* no top clip -> full region */
  if(y==h) {
    *t = *b = 0;
    fully = 1;
  } else {
    fully = 0;
    /* find bot clip */
    for(y=yoffset;y<h-yoffset;y++) {
      if(yd[h-y-1]>=ythres) {
	*b = y;
	break;
      }
    }
    /* same threshold -> no clipping useful */
    if(*b + *t == height) {
      *t = *b = 0;
      fully = 1;
    }
  }

  /* ----- horizontal operation ----- */
  /* find left clip */
  for(x=xoffset;x<w-xoffset;x++) {
    if(xd[x]>=xthres) {
      *l = x;
      break;
    }
  }
  /* no left clip -> use full region */
  if(x==w) {
    *l = *r = 0;
    fullx = 1;
  } else {
    fullx = 0;
    /* find bot clip */
    for(x=xoffset;x<w-xoffset;x++) {
      if(xd[w-x-1]>=xthres) {
	*r = x;
	break;
      }
    }
    /* same threshold -> no clipping useful */
    if(*l + *r == width) {
      *l = *r = 0;
      fullx = 1;
    }
  }

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
      saveGnuplot(file,xd,w,*l,w-1-*r,xoffset,xthres);
      file[len+1] = 'y';
      saveGnuplot(file,yd,h,*t,h-1-*b,yoffset,ythres);
      free(file);
    }
  }

  free(xd);
  free(yd);

  return(!(fullx && fully));
}

/* read the next header line of pgm file */
unsigned char *readHeaderLine(FILE *fh)
{
  static unsigned char buf[1024];

  do {
    /* read error */
    if(fgets(buf,1023,fh)==0) {
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
  if(strcmp("P5\n",line)!=0) {
    fprintf(stderr,"No P5 pgm!\n");
    return 0;
  }
  /* read image dimension */
  line = readHeaderLine(fh);
  if(sscanf(line,"%d %d %d",width,height,&dummy)!=3) {
    if(sscanf(line,"%d %d",width,height)!=2) {
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
  ptr1 = memory +  t    * w + l;
  ptr2 = memory + (t+1) * w - 1 - r;
  for(y=t;y<(h-b);y++) {
    *ptr1 = 255 - *ptr1;
    *ptr2 = 255 - *ptr2;
    ptr1 += w;
    ptr2 += w;
  }

  /* draw horizontal markers */
  ptr1 = memory +  t      * w + l + 1;
  ptr2 = memory + (h-1-b) * w + l + 1;
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
" -t <thres>[,<ythres>]     set threshold values    (default: %d,%d)\n"
" -s <safety>[,<ysafety>]   add a safety border     (default: none)\n"
" -b <modulo>[,<ymodulo>]   align clip borders      (default: none)\n"
" -f <modulo>[,<ymodulo>]   align output frame size (default: none)\n"
" -o <offset>[,<yoffset>]   search begin offset     (default: 0 0)\n"
" -e                        expand and do not shrink in frame alignment\n"
" -y                        input pgm files are yuv files from mplayer\n"
" -v                        enable verbose mode\n"
" -p                        plot results to *x.eps, *y.eps (req. gnuplot)\n"
" -w                        write PGM files with border markers (*-m.pgm)\n"
	  ,xthres,ythres);
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
  int top,bottom,left,right;
  int width,height;
  int count;

  /* parse args */
  int c;
  while((c=getopt(argc,argv,"t:s:f:b:o:veypw"))!=-1) {
    switch(c) {
    case 't':
      twoValArg(optarg,&xthres,&ythres);
      break;
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

#define MAX 2048

  top    = MAX;
  bottom = MAX;
  left   = MAX;
  right  = MAX;

  /* loop over images */
  count = 0;
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

      /* keep largest region */
      if(l<left)
	left = l;
      if(r<right)
	right = r;
      if(t<top)
	top = t;
      if(b<bottom)
	bottom = b;

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
  
  /* consistency check */
  if(top > (height-1-bottom)) {
    fprintf(stderr,"Inconsistency along Y - ignoring values!\n");
    top = bottom = 0;
  }
  if(left > (width-1-right)) {
    fprintf(stderr,"Inconsistency along X - ignoring values!\n");
    left = right = 0;
  }

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

  /* generate transcode-friendly result */
  printf("%d,%d,%d,%d\n",top,left,bottom,right);

  exit(0);
}
