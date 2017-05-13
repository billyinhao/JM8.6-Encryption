/**********************************************************************
 * Software Copyright Licensing Disclaimer
 *
 * This software module was originally developed by contributors to the
 * course of the development of ISO/IEC 14496-10 for reference purposes
 * and its performance may not have been optimized.  This software
 * module is an implementation of one or more tools as specified by
 * ISO/IEC 14496-10.  ISO/IEC gives users free license to this software
 * module or modifications thereof. Those intending to use this software
 * module in products are advised that its use may infringe existing
 * patents.  ISO/IEC have no liability for use of this software module
 * or modifications thereof.  The original contributors retain full
 * rights to modify and use the code for their own purposes, and to
 * assign or donate the code to third-parties.
 *
 * This copyright notice must be included in all copies or derivative
 * works.  Copyright (c) ISO/IEC 2004.
 **********************************************************************/

/*!
 ***********************************************************************
 * \file macroblock.c
 *
 * \brief
 *     Decode a Macroblock
 *
 * \author
 *    Main contributors (see contributors.h for copyright, address and affiliation details)
 *    - Inge Lille-Lang鴜               <inge.lille-langoy@telenor.com>
 *    - Rickard Sjoberg                 <rickard.sjoberg@era.ericsson.se>
 *    - Jani Lainema                    <jani.lainema@nokia.com>
 *    - Sebastian Purreiter             <sebastian.purreiter@mch.siemens.de>
 *    - Thomas Wedi                     <wedi@tnt.uni-hannover.de>
 *    - Detlev Marpe                    <marpe@hhi.de>
 *    - Gabi Blaettermann               <blaetter@hhi.de>
 *    - Ye-Kui Wang                     <wyk@ieee.org>
 ***********************************************************************
*/

#include "contributors.h"

#include <math.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

#include "global.h"
#include "mbuffer.h"
#include "elements.h"
#include "errorconcealment.h"
#include "macroblock.h"
#include "fmo.h"
#include "cabac.h"
#include "vlc.h"
#include "image.h"
#include "mb_access.h"
#include "biaridecod.h"

#if TRACE
#define TRACE_STRING(s) strncpy(currSE.tracestring, s, TRACESTRING_SIZE)
#else
#define TRACE_STRING(s) // do nothing
#endif

extern int last_dquant;
extern ColocatedParams *Co_located;
//////////////////////////////////////////////////////////////////////////////////////////////////
extern double u, init_con;
extern int cn;
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


static void SetMotionVectorPredictor (struct img_par  *img,
                                      int             *pmv_x,
                                      int             *pmv_y,
                                      int             ref_frame,
                                      int             list,
                                      int             ***refPic,
                                      int             ****tmp_mv,
                                      int             block_x,
                                      int             block_y,
                                      int             blockshape_x,
                                      int             blockshape_y);


/*!
 ************************************************************************
 * \brief
 *    initializes the current macroblock
 ************************************************************************
 */
 //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

short logistic_generate(short mod)
{
	double temp;
	 init_con = u*init_con*(1-init_con);
	       if(init_con == 0.5) 
		     init_con = 0.4;
		   temp=init_con;
		   
   return ((long long)(temp*1.0e+14)) % mod;
}
static int value_encryption(int  value,int mod)
{
/////////////////////////////////////////
  //return value^logistic_generate(mod);
  //////////////////////////////////////
  cn=value^logistic_generate(mod);
  return cn;
}
 void show_matrix_value(struct img_par *img)
{
  int i,j;
  putchar('\n');
   for(j=0;j<4;j++)
   {
	  for(i=0;i<4;i++)
	  {
	   printf("%d ",img->ipredmode[img->block_x+i][img->block_y+j]);
	  }
	  putchar('\n');
   }
}

void decryption_intra_mode(struct img_par *img)
{
   static int count=0;
   int i,j;
   int location_of_macroblock=0;
   if(img->current_mb_nr==0)
   location_of_macroblock=1;
   else if(img->current_mb_nr<img->PicWidthInMbs)
   location_of_macroblock=2;
   else if(img->current_mb_nr%img->PicWidthInMbs==0)
   location_of_macroblock=3;
   //if(count++<1)
	if(location_of_macroblock==0)
	{
	  for(j=0;j<4;j++)
	   for(i=0;i<4;i++)
	   if(img->ipredmode[img->block_x+i][img->block_y+j]!=-1)
		img->ipredmode[img->block_x+i][img->block_y+j]=value_encryption(img->ipredmode[img->block_x+i][img->block_y+j],8);
	}
	 else if(location_of_macroblock==1)
	{
	 //show_matrix_value(img);
	 for(j=1;j<4;j++)
	  for(i=1;i<4;i++)
	  {
	   if(img->ipredmode[img->block_x+i][img->block_y+j]!=-1)
		img->ipredmode[img->block_x+i][img->block_y+j]=value_encryption(img->ipredmode[img->block_x+i][img->block_y+j],8);
	  }
	}
	 else if(location_of_macroblock==2)
	{
	 for(j=1;j<4;j++)
	  for(i=0;i<4;i++)
	   if(img->ipredmode[img->block_x+i][img->block_y+j]!=-1)
		  img->ipredmode[img->block_x+i][img->block_y+j]=value_encryption(img->ipredmode[img->block_x+i][img->block_y+j],8);
	}
	 else if(location_of_macroblock==3)
	{
	 for(j=0;j<4;j++)
	  for(i=1;i<4;i++)
	   if(img->ipredmode[img->block_x+i][img->block_y+j]!=-1)
		  img->ipredmode[img->block_x+i][img->block_y+j]=value_encryption(img->ipredmode[img->block_x+i][img->block_y+j],8);
	}
}
int decryption_C_intra_mode(int  c_ipred)
{
   int location_of_macroblock=0;
   if(img->current_mb_nr==0)
   location_of_macroblock=1;
   else if(img->current_mb_nr<img->PicWidthInMbs)
   location_of_macroblock=2;
   else if(img->current_mb_nr%img->PicWidthInMbs==0)
   location_of_macroblock=3;
    return (location_of_macroblock==0)? value_encryption(c_ipred,4):c_ipred;
}
static void DCT_no_coeff_decryption(int * pLevel,int * Run,int cdc,int DCT)
{
  int level,p;
  int numtrailingones,lastcoeff;
  int numcoeff,shift;
  int vlcnum,level_two_or_higher;
  int levabs;
  int levelcode;
  const int incVlc[] = {0,3,6,12,24,48,32768}; 
  int numPrefix,numSuffix;
  int escape;
  level=1;
  numtrailingones=0;
  lastcoeff = 0;
  numcoeff =0;
  for(p=0;(p<(cdc?4:16))&&level!=0;p++)
  {
    level = pLevel[p]; 
	if(level)
	{
      if(abs(level)==1)
	  {
	   numtrailingones++;
	   if(numtrailingones>3)
		numtrailingones=3;
	  }
	  else
	   numtrailingones=0;
	   numcoeff++;
	   lastcoeff=p;
    }
  }
  level_two_or_higher = 1;//标识拖尾系数的值，如果为1代表拖尾系数小于3，第一个Level的计算要加1(>0)或减1(<0)
  if (numcoeff > 3 && numtrailingones == 3)
   level_two_or_higher = 0;
  if (numcoeff > 10 && numtrailingones < 3)
   vlcnum = 1;//suffixLength????
  else
   vlcnum = 0;
  //拖尾系数处理
  for(p=lastcoeff;p>lastcoeff - numtrailingones;p--)
 {
  if(logistic_generate(2)==0)
   pLevel[p]*=-1;
 }
  //对非拖尾系数进行处理
 if(DCT>1)
 for (p = lastcoeff - numtrailingones; p >= 0; p--)
 {
  level=pLevel[p];
  //只处理vlcnum(suffixlength)>0
  if(vlcnum!=0)
  {
   shift=vlcnum-1;
   escape=(15<<shift)+1;
   levabs=abs(level);
   if(levabs<escape)
   {
	levelcode=level<0?2*levabs-1:2*levabs-2;
	numPrefix=levelcode/(1<<vlcnum);
	if(DCT>2)
	{
	if(numPrefix<=2&&vlcnum<=5)
	{
	 if(p==lastcoeff - numtrailingones&&level_two_or_higher==1)
      numPrefix=(numPrefix-1+2-logistic_generate(2))%2+1;
	 else
	  numPrefix=(numPrefix+3-logistic_generate(3))%3;
	}
	if(DCT>3)
	{
	if(numPrefix>2&&vlcnum<=5)
	{
	  numPrefix=(numPrefix-3+12-logistic_generate(12))%12+3;
	}
	if(vlcnum==6)
	{
	  numPrefix=value_encryption(numPrefix,15);
	}
	}
	}
    numSuffix=levelcode%(1<<vlcnum);
	numSuffix=value_encryption(numSuffix,1<<vlcnum);
	levelcode=numPrefix*(1<<vlcnum)+numSuffix;
	pLevel[p]=(levelcode%2==0)?(levelcode+2)>>1:-((levelcode+1)>>1);
   }
  }
	if (abs(pLevel[p])>incVlc[vlcnum])
     vlcnum++;
    if (p == lastcoeff - numtrailingones && abs(pLevel[p])>3)
     vlcnum = 2;
 }

}
int MVD_decryption(int mvd,int cavlc)
{ 
  int absmvd,i,j,sub=0,temp=1,codenum;
  absmvd=abs(mvd);
  codenum=(mvd>0)?2*absmvd-1:2*absmvd;
  if(cavlc)
  {
    absmvd=abs(mvd);
	for(i=0;i<16&&absmvd!=0;i++)
	{
		absmvd/=2;
	}
	for(j=0;j<=i-1;j++)
	{
		sub+=temp;
		temp*=2;
	}
	if(i!=0)
	{
	 codenum=value_encryption(codenum-sub,temp);
	 codenum=sub+codenum;
	}
	return (!(codenum%2))?-(codenum/2):(codenum+1)/2;
  }
  else
  return -mvd;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void start_macroblock(struct img_par *img,struct inp_par *inp, int CurrentMBInScanOrder)
{
  int i,j,k,l;
  Macroblock *currMB;   // intialization code deleted, see below, StW

  assert (img->current_mb_nr >=0 && img->current_mb_nr < img->PicSizeInMbs);

  currMB = &img->mb_data[img->current_mb_nr];

  /* Update coordinates of the current macroblock */
  if (img->MbaffFrameFlag)
  {
    img->mb_x = (img->current_mb_nr)%((2*img->width)/MB_BLOCK_SIZE);
    img->mb_y = 2*((img->current_mb_nr)/((2*img->width)/MB_BLOCK_SIZE));

    if (img->mb_x % 2)
    {
      img->mb_y++;
    }

    img->mb_x /= 2;
  }
  else
  {
    img->mb_x = (img->current_mb_nr)%(img->width/MB_BLOCK_SIZE);
    img->mb_y = (img->current_mb_nr)/(img->width/MB_BLOCK_SIZE);
  }
  
  /* Define vertical positions */
  img->block_y = img->mb_y * BLOCK_SIZE;      /* luma block position */
  img->pix_y   = img->mb_y * MB_BLOCK_SIZE;   /* luma macroblock position */
  img->pix_c_y = img->mb_y * MB_BLOCK_SIZE/2; /* chroma macroblock position */
  
  /* Define horizontal positions */
  img->block_x = img->mb_x * BLOCK_SIZE;      /* luma block position */
  img->pix_x   = img->mb_x * MB_BLOCK_SIZE;   /* luma pixel position */
  img->pix_c_x = img->mb_x * MB_BLOCK_SIZE/2; /* chroma pixel position */

  // Save the slice number of this macroblock. When the macroblock below
  // is coded it will use this to decide if prediction for above is possible
  currMB->slice_nr = img->current_slice_nr;

  if (img->current_slice_nr >= MAX_NUM_SLICES)
  {
    error ("maximum number of supported slices exceeded, please recompile with increased value for MAX_NUM_SLICES", 200);
  }

  dec_picture->slice_id[img->mb_x][img->mb_y] = img->current_slice_nr;
  if (img->current_slice_nr > dec_picture->max_slice_id)
  {
    dec_picture->max_slice_id=img->current_slice_nr;
  }
  
  CheckAvailabilityOfNeighbors(img);

  // Reset syntax element entries in MB struct
  currMB->qp          = img->qp ;
  currMB->mb_type     = 0;
  currMB->delta_quant = 0;
  currMB->cbp         = 0;
  currMB->cbp_blk     = 0;
  currMB->c_ipred_mode= DC_PRED_8; //GB

  for (l=0; l < 2; l++)
    for (j=0; j < BLOCK_MULTIPLE; j++)
      for (i=0; i < BLOCK_MULTIPLE; i++)
        for (k=0; k < 2; k++)
          currMB->mvd[l][j][i][k] = 0;

  currMB->cbp_bits   = 0;

  // initialize img->m7 for ABT
  for (j=0; j<MB_BLOCK_SIZE; j++)
    for (i=0; i<MB_BLOCK_SIZE; i++)
      img->m7[i][j] = 0;

  // store filtering parameters for this MB 
  currMB->LFDisableIdc = img->currentSlice->LFDisableIdc;
  currMB->LFAlphaC0Offset = img->currentSlice->LFAlphaC0Offset;
  currMB->LFBetaOffset = img->currentSlice->LFBetaOffset;

}

/*!
 ************************************************************************
 * \brief
 *    set coordinates of the next macroblock
 *    check end_of_slice condition 
 ************************************************************************
 */
int exit_macroblock(struct img_par *img,struct inp_par *inp,int eos_bit)
{

  //! The if() statement below resembles the original code, which tested 
  //! img->current_mb_nr == img->PicSizeInMbs.  Both is, of course, nonsense
  //! In an error prone environment, one can only be sure to have a new
  //! picture by checking the tr of the next slice header!

// printf ("exit_macroblock: FmoGetLastMBOfPicture %d, img->current_mb_nr %d\n", FmoGetLastMBOfPicture(), img->current_mb_nr);
  img->num_dec_mb++;
  
  if (img->num_dec_mb == img->PicSizeInMbs)
//  if (img->current_mb_nr == FmoGetLastMBOfPicture(currSlice->structure))
  {
//thb
/*
    if (currSlice->next_header != EOS)
      currSlice->next_header = SOP;
*/
//the
    assert (nal_startcode_follows (img, inp, eos_bit) == TRUE);
    return TRUE;
  }
  // ask for last mb in the slice  UVLC
  else
  {
// printf ("exit_macroblock: Slice %d old MB %d, now using MB %d\n", img->current_slice_nr, img->current_mb_nr, FmoGetNextMBNr (img->current_mb_nr));

    img->current_mb_nr = FmoGetNextMBNr (img->current_mb_nr);

    if (img->current_mb_nr == -1)     // End of Slice group, MUST be end of slice
    {
      assert (nal_startcode_follows (img, inp, eos_bit) == TRUE);
      return TRUE;
    }

    if(nal_startcode_follows(img, inp, eos_bit) == FALSE) 
      return FALSE;

    if(img->type == I_SLICE  || img->type == SI_SLICE || active_pps->entropy_coding_mode_flag == CABAC)
      return TRUE;
    if(img->cod_counter<=0)
      return TRUE;
    return FALSE;
  }
}

/*!
 ************************************************************************
 * \brief
 *    Interpret the mb mode for P-Frames
 ************************************************************************
 */
void interpret_mb_mode_P(struct img_par *img)
{
  int i;
  const int ICBPTAB[6] = {0,16,32,15,31,47};	//++ CBP的6个比特位的十进制表示
  Macroblock *currMB = &img->mb_data[img->current_mb_nr];
  int         mbmode = currMB->mb_type;

#define ZERO_P8x8     (mbmode==5)
#define MODE_IS_P8x8  (mbmode==4 || mbmode==5)
#define MODE_IS_I4x4  (mbmode==6)
#define I16OFFSET     (mbmode-7)
#define MODE_IS_IPCM  (mbmode==31)

  if(mbmode <4)	//++ P_Skip、帧间宏块分割
  {
    currMB->mb_type = mbmode;
    for (i=0;i<4;i++)
    {
      currMB->b8mode[i]   = mbmode;
      currMB->b8pdir[i]   = 0;
    }
  }
  else if(MODE_IS_P8x8)	//++ 帧间子宏块分割
  {
    currMB->mb_type = P8x8;
    img->allrefzero = ZERO_P8x8;
  }
  else if(MODE_IS_I4x4)	//++ 帧内4*4
  {
    currMB->mb_type = I4MB;
    for (i=0;i<4;i++)
    {
      currMB->b8mode[i] = IBLOCK;
      currMB->b8pdir[i] = -1;
    }
  }
  else if(MODE_IS_IPCM)	//++ I_PCM
  {
    currMB->mb_type=IPCM;
    
    for (i=0;i<4;i++) 
    {
      currMB->b8mode[i]=0; currMB->b8pdir[i]=-1; 
    }
    currMB->cbp= -1;
    currMB->i16mode = 0;
  }
  else	//++ 帧内16*16
  {
    currMB->mb_type = I16MB;
    for (i=0;i<4;i++) {currMB->b8mode[i]=0; currMB->b8pdir[i]=-1; }
    currMB->cbp= ICBPTAB[(I16OFFSET)>>2];
    currMB->i16mode = (I16OFFSET) & 0x03;
  }
}

/*!
 ************************************************************************
 * \brief
 *    Interpret the mb mode for I-Frames
 ************************************************************************
 */
void interpret_mb_mode_I(struct img_par *img)
{
  int i;
  const int ICBPTAB[6] = {0,16,32,15,31,47};	//++ CBP的6个比特位的十进制表示
  Macroblock *currMB   = &img->mb_data[img->current_mb_nr];
  int         mbmode   = currMB->mb_type;

  if (mbmode==0)	//++ 帧内4*4
  {
    currMB->mb_type = I4MB;
    for (i=0;i<4;i++) {currMB->b8mode[i]=IBLOCK; currMB->b8pdir[i]=-1; }
  }
  else if(mbmode==25)	//++ I_PCM
  {
    for (i=0;i<4;i++) {currMB->b8mode[i]=0; currMB->b8pdir[i]=-1; }
    currMB->cbp= -1;
    currMB->i16mode = 0;

  }
  else	//++ 帧内16*16
  {
    currMB->mb_type = I16MB;
    for (i=0;i<4;i++) {currMB->b8mode[i]=0; currMB->b8pdir[i]=-1; }
    currMB->cbp= ICBPTAB[(mbmode-1)>>2];
    currMB->i16mode = (mbmode-1) & 0x03;
  }
}

/*!
 ************************************************************************
 * \brief
 *    Interpret the mb mode for B-Frames
 ************************************************************************
 */
void interpret_mb_mode_B(struct img_par *img)
{
  static const int offset2pdir16x16[12]   = {0, 0, 1, 2, 0,0,0,0,0,0,0,0};
  static const int offset2pdir16x8[22][2] = {{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{1,1},{0,0},{0,1},{0,0},{1,0},
                                             {0,0},{0,2},{0,0},{1,2},{0,0},{2,0},{0,0},{2,1},{0,0},{2,2},{0,0}};
  static const int offset2pdir8x16[22][2] = {{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{1,1},{0,0},{0,1},{0,0},
                                             {1,0},{0,0},{0,2},{0,0},{1,2},{0,0},{2,0},{0,0},{2,1},{0,0},{2,2}};

  const int ICBPTAB[6] = {0,16,32,15,31,47};
  Macroblock *currMB = &img->mb_data[img->current_mb_nr];

  int i, mbmode;
  int mbtype  = currMB->mb_type;
  int *b8mode = currMB->b8mode;
  int *b8pdir = currMB->b8pdir;

  if (mbtype == 3 || (mbtype>15 && mbtype <22))
  {
	  mbtype = mbtype;
  }
  //--- set mbtype, b8type, and b8pdir ---
  if (mbtype==0)       //++ B_Skip、B_Direct_16*16
  {
      mbmode=0;       for(i=0;i<4;i++) {b8mode[i]=0;          b8pdir[i]=2; }
  }
  else if (mbtype==23) //++ intra4x4
  {
    mbmode=I4MB;    for(i=0;i<4;i++) {b8mode[i]=IBLOCK;     b8pdir[i]=-1; }
  }
  else if ((mbtype>23) && (mbtype<48) ) // intra16x16
  {
    mbmode=I16MB;   for(i=0;i<4;i++) {b8mode[i]=0;          b8pdir[i]=-1; }
    currMB->cbp     = ICBPTAB[(mbtype-24)>>2];
    currMB->i16mode = (mbtype-24) & 0x03;
  }
  else if (mbtype==22) //++ inter8x8(+split)
  {
    mbmode=P8x8;       // b8mode and pdir is transmitted in additional codewords
  }
  else if (mbtype<4)   //++ inter16x16
  {
    mbmode=1;       for(i=0;i<4;i++) {b8mode[i]=1;          b8pdir[i]=offset2pdir16x16[mbtype]; }
  }
  else if(mbtype==48)	//++ I_PCM
  {
    mbmode=IPCM;
    for (i=0;i<4;i++) {currMB->b8mode[i]=0; currMB->b8pdir[i]=-1; }
    currMB->cbp= -1;
    currMB->i16mode = 0;
  }

  else if (mbtype%2==0) //++ inter16x8
  {
    mbmode=2;       for(i=0;i<4;i++) {b8mode[i]=2;          b8pdir[i]=offset2pdir16x8 [mbtype][i/2]; }
  }
  else	//++ inter8*16
  {
    mbmode=3;       for(i=0;i<4;i++) {b8mode[i]=3;          b8pdir[i]=offset2pdir8x16 [mbtype][i%2]; }
  }
  currMB->mb_type = mbmode;
}
/*!
 ************************************************************************
 * \brief
 *    Interpret the mb mode for SI-Frames
 ************************************************************************
 */
void interpret_mb_mode_SI(struct img_par *img)
{
  int i;
  const int ICBPTAB[6] = {0,16,32,15,31,47};
  Macroblock *currMB   = &img->mb_data[img->current_mb_nr];
  int         mbmode   = currMB->mb_type;

  if (mbmode==0)
  {
    currMB->mb_type = SI4MB;
    for (i=0;i<4;i++) {currMB->b8mode[i]=IBLOCK; currMB->b8pdir[i]=-1; }
    img->siblock[img->mb_x][img->mb_y]=1;
  }
  else if (mbmode==1)
  {
    currMB->mb_type = I4MB;
    for (i=0;i<4;i++) {currMB->b8mode[i]=IBLOCK; currMB->b8pdir[i]=-1; }
  }
  else if(mbmode==26)
  {
    currMB->mb_type=IPCM;

    for (i=0;i<4;i++) {currMB->b8mode[i]=0; currMB->b8pdir[i]=-1; }
    currMB->cbp= -1;
    currMB->i16mode = 0;

  }

  else
  {
    currMB->mb_type = I16MB;
    for (i=0;i<4;i++) {currMB->b8mode[i]=0; currMB->b8pdir[i]=-1; }
    currMB->cbp= ICBPTAB[(mbmode-1)>>2];
    currMB->i16mode = (mbmode-2) & 0x03;
  }
}
/*!
 ************************************************************************
 * \brief
 *    init macroblock I and P frames
 ************************************************************************
 */
void init_macroblock(struct img_par *img)
{
  int i,j;

  for (i=0;i<BLOCK_SIZE;i++)
  {                           // reset vectors and pred. modes
    for(j=0;j<BLOCK_SIZE;j++)
    {
      dec_picture->mv[LIST_0][img->block_x+i][img->block_y+j][0]=0;
      dec_picture->mv[LIST_0][img->block_x+i][img->block_y+j][1]=0;
      dec_picture->mv[LIST_1][img->block_x+i][img->block_y+j][0]=0;
      dec_picture->mv[LIST_1][img->block_x+i][img->block_y+j][1]=0;

      img->ipredmode[img->block_x+i][img->block_y+j] = DC_PRED;
    }
  }

  for (j=0; j<BLOCK_SIZE; j++)
    for (i=0; i<BLOCK_SIZE; i++)
    {
      dec_picture->ref_idx[LIST_0][img->block_x+i][img->block_y+j] = -1;
      dec_picture->ref_idx[LIST_1][img->block_x+i][img->block_y+j] = -1;
      dec_picture->ref_pic_id[LIST_0][img->block_x+i][img->block_y+j] = INT64_MIN;
      dec_picture->ref_pic_id[LIST_1][img->block_x+i][img->block_y+j] = INT64_MIN;
    }
}


/*!
 ************************************************************************
 * \brief
 *    Sets mode for 8x8 block
 ************************************************************************
 */
void SetB8Mode (struct img_par* img, Macroblock* currMB, int value, int i)
{
  static const int p_v2b8 [ 5] = {4, 5, 6, 7, IBLOCK};
  static const int p_v2pd [ 5] = {0, 0, 0, 0, -1};
  static const int b_v2b8 [14] = {0, 4, 4, 4, 5, 6, 5, 6, 5, 6, 7, 7, 7, IBLOCK};	//++ B片子宏块类型。0：Direct_8*8，4：8*8，5：8*4，6：4*8，7：4*4，IBLOCK：帧内
  static const int b_v2pd [14] = {2, 0, 1, 2, 0, 0, 1, 1, 2, 2, 0, 1, 2, -1};	//++ B片子宏块分割的预测方向。0：前向，1：后向，2：双向（或Direct），-1：无预测（帧内）

  if (img->type==B_SLICE)
  {
    currMB->b8mode[i]   = b_v2b8[value];
    currMB->b8pdir[i]   = b_v2pd[value];

  }
  else
  {
    currMB->b8mode[i]   = p_v2b8[value];
    currMB->b8pdir[i]   = p_v2pd[value];
  }

}


void reset_coeffs()
{
  int i, j, iii, jjj;

  // reset luma coeffs
  for (i=0;i<BLOCK_SIZE;i++)
  { 
    for (j=0;j<BLOCK_SIZE;j++)
      for(iii=0;iii<BLOCK_SIZE;iii++)
        for(jjj=0;jjj<BLOCK_SIZE;jjj++)
          img->cof[i][j][iii][jjj]=0;
  }

  // reset chroma coeffs
  for (j=4;j<6;j++)
  { 
    for (i=0;i<4;i++)
      for (iii=0;iii<4;iii++)
        for (jjj=0;jjj<4;jjj++)
          img->cof[i][j][iii][jjj]=0;
  }
  
  // CAVLC
  for (i=0; i < 4; i++)
    for (j=0; j < 6; j++)
      img->nz_coeff[img->current_mb_nr][i][j]=0;  

}

void field_flag_inference()
{
  Macroblock *currMB = &img->mb_data[img->current_mb_nr];

  if (currMB->mbAvailA)
  {
    currMB->mb_field = img->mb_data[currMB->mbAddrA].mb_field;
  }
  else
  {
    // check top macroblock pair
    if (currMB->mbAvailB)
    {
      currMB->mb_field = img->mb_data[currMB->mbAddrB].mb_field;
    }
    else
      currMB->mb_field = 0;
  }
  
}

/*!
 ************************************************************************
 * \brief
 *    Get the syntax elements from the NAL
 ************************************************************************
 */

/************************************************************************************************************
read_one_macroblock 函数与 7.3.4 有一定对应关系,看看 7.3.4 中粗体字 mb_skip_run（mb_skip_run 表示从当前
宏块开始连续出现的 skip 宏块的数量）,可以看见该语法元素存在的条件是if( !entropy_coding_mode_flag )也就是采
用 CAVLC 编码时候,这里只以 CAVLC 为例进行讲解。解码器在这里解码出 mb_skip_run 后，就会根据 mb_skip_run 的值
进行相应处理（mb_skip_run 在代码中是变量 img->cod_counter ）。img->cod_counter 大于 0 就直接调用处理 skip 宏
块的代码 
*****************************************************************************************************************/ 

/****************************************************************************************************************
read_one_macroblocko函数首先初始化一些参数和变量。接着根据宏块类型的不同采用不同的interpret_mb_mode()函数来对宏块
进行解析，如I、P、B、sI宏块。然后判断该宏块是否为8x8子宏块分割，如果是，则解析子宏块分割模式。init_macroblock()函
数主要作用是对亮度块中16个4x4块的预测模式、运动矢量进行初始化。然后IS_DIRECT判断该宏块是否B_kip类型宏块，如果是，
则按相应算法处理。接着IS_COPY判断该宏块是否COPY类型宏块，如果是，则按相应算法处理。下一步，read_iPredes_modes从比
特流中读取亮度块和色度块的预测模式。以亮度块为例，读取时是以4x4块为单位读取的（有一个读取顺序，见标准）。
IS_INTERMV 判断该宏块是否帧间编码宏块，如果是，则通过函数readMotionInfoFromNAl()从码流中读取运动矢量残差。最后是分
别对亮度和色度宏块进行CAVLC嫡解码。解码最小单位是4x4块，每个4x4块在16x16宏块中解码顺序见标准。
****************************************************************************************************************/
 
int read_one_macroblock(struct img_par *img,struct inp_par *inp)
{
  int i;

  SyntaxElement currSE;
  Macroblock *currMB = &img->mb_data[img->current_mb_nr];

  Slice *currSlice = img->currentSlice;
  DataPartition *dP;
  int *partMap = assignSE2partition[currSlice->dp_mode];
  Macroblock *topMB = NULL;
  int  prevMbSkipped = 0;
  int  img_block_y;
  int  check_bottom, read_bottom, read_top;
  
  if (img->MbaffFrameFlag)
  {
    if (img->current_mb_nr%2)
    {
      topMB= &img->mb_data[img->current_mb_nr-1];
      if(!(img->type == B_SLICE))
        prevMbSkipped = (topMB->mb_type == 0);
      else 
        prevMbSkipped = (topMB->mb_type == 0 && topMB->cbp == 0);
    }
    else 
      prevMbSkipped = 0;
  }
  
  if (img->current_mb_nr%2 == 0)
    currMB->mb_field = 0;
  else
    currMB->mb_field = img->mb_data[img->current_mb_nr-1].mb_field;


  currMB->qp = img->qp ;

  currSE.type = SE_MBTYPE;

  //  read MB mode *****************************************************************
  dP = &(currSlice->partArr[partMap[currSE.type]]);	//++ 对currSlice->partArr分配内存空间在decod.c文件malloc_slice()函数最后一行
													//++ 对其赋值是在image.c文件read_new_slice()函数while循环的switch语句中
  if (active_pps->entropy_coding_mode_flag == UVLC || dP->bitstream->ei_flag)   currSE.mapping = linfo_ue;

  if(img->type == I_SLICE || img->type == SI_SLICE)
  {
    // read MB aff
    if (img->MbaffFrameFlag && img->current_mb_nr%2==0)
    {
      TRACE_STRING("mb_field_decoding_flag");
      if (active_pps->entropy_coding_mode_flag == UVLC || dP->bitstream->ei_flag)
      {
        currSE.len = 1;
        readSyntaxElement_FLC(&currSE, dP->bitstream);
      }
      else
      {
        currSE.reading = readFieldModeInfo_CABAC;
        dP->readSyntaxElement(&currSE,img,inp,dP);
      }
      currMB->mb_field = currSE.value1;
    }
    if(active_pps->entropy_coding_mode_flag  == CABAC)
      CheckAvailabilityOfNeighborsCABAC();

    //  read MB type
    TRACE_STRING("mb_type");
    currSE.reading = readMB_typeInfo_CABAC;
    dP->readSyntaxElement(&currSE,img,inp,dP);

    currMB->mb_type = currSE.value1;
    if(!dP->bitstream->ei_flag)
      currMB->ei_flag = 0;
  }

  // non I/SI-slice CABAC
  else if (active_pps->entropy_coding_mode_flag == CABAC)
  {
    // read MB skipflag,从码流中读取mb_skip_run这个语法元素 
    if (img->MbaffFrameFlag && (img->current_mb_nr%2 == 0||prevMbSkipped))
      field_flag_inference();
    
    CheckAvailabilityOfNeighborsCABAC();
    TRACE_STRING("mb_skip_flag");
    currSE.reading = readMB_skip_flagInfo_CABAC;
    dP->readSyntaxElement(&currSE,img,inp,dP);
    currMB->mb_type = currSE.value1;

    if (img->type==B_SLICE)
      currMB->cbp = currSE.value2;
    if(!dP->bitstream->ei_flag)
      currMB->ei_flag = 0;
    
    if ((img->type==B_SLICE) && currSE.value1==0 && currSE.value2==0)
      img->cod_counter=0;
      
    // read MB aff
      if (img->MbaffFrameFlag) 
      {
        check_bottom=read_bottom=read_top=0;
        if (img->current_mb_nr%2==0)
        {
          check_bottom =  (img->type!=B_SLICE)? 
            (currMB->mb_type == 0):
          (currMB->mb_type == 0 && currMB->cbp == 0);
          read_top = !check_bottom;
        }
        else
          read_bottom = (img->type!=B_SLICE)? 
          (topMB->mb_type == 0 && currMB->mb_type != 0) :
        ((topMB->mb_type == 0 && topMB->cbp == 0) && (currMB->mb_type != 0 || currMB->cbp != 0));
        
        if (read_bottom || read_top)
        {
          TRACE_STRING("mb_field_decoding_flag");
          currSE.reading = readFieldModeInfo_CABAC;
          dP->readSyntaxElement(&currSE,img,inp,dP);
          currMB->mb_field = currSE.value1;
        }
        if (check_bottom)
          check_next_mb_and_get_field_mode_CABAC(&currSE,img,inp,dP);
        
      }
      if(active_pps->entropy_coding_mode_flag  == CABAC)
        CheckAvailabilityOfNeighborsCABAC();
      
    // read MB type
    if (currMB->mb_type != 0 )
    {
      currSE.reading = readMB_typeInfo_CABAC;
      TRACE_STRING("mb_type");
      dP->readSyntaxElement(&currSE,img,inp,dP);
      currMB->mb_type = currSE.value1;
      if(!dP->bitstream->ei_flag)
        currMB->ei_flag = 0;
    }
  }
  // VLC Non-Intra
  else
  {
    if(img->cod_counter == -1)
    {
      TRACE_STRING("mb_skip_run");
      dP->readSyntaxElement(&currSE,img,inp,dP);
      img->cod_counter = currSE.value1;
    }
    if (img->cod_counter==0)
    {
      // read MB aff
      if ((img->MbaffFrameFlag) && ((img->current_mb_nr%2==0) || ((img->current_mb_nr%2) && prevMbSkipped)))
      {
        TRACE_STRING("mb_field_decoding_flag");
        currSE.len = 1;
        readSyntaxElement_FLC(&currSE, dP->bitstream);
        currMB->mb_field = currSE.value1;
      }
      
      // read MB type
      TRACE_STRING("mb_type");
      dP->readSyntaxElement(&currSE,img,inp,dP);
      if(img->type == P_SLICE || img->type == SP_SLICE)
        currSE.value1++;
      currMB->mb_type = currSE.value1;
      if(!dP->bitstream->ei_flag)
        currMB->ei_flag = 0;
      img->cod_counter--;
    } 
    else
    {
      img->cod_counter--;
      currMB->mb_type = 0;
      currMB->ei_flag = 0;

      // read field flag of bottom block
      if(img->MbaffFrameFlag)
      {
        if(img->cod_counter == 0 && (img->current_mb_nr%2 == 0))
        {
          TRACE_STRING("mb_field_decoding_flag (of coded bottom mb)");
          currSE.len = 1;
          readSyntaxElement_FLC(&currSE, dP->bitstream);
          dP->bitstream->frame_bitoffset--;
          currMB->mb_field = currSE.value1;
        }
        else if(img->cod_counter > 0 && (img->current_mb_nr%2 == 0))
        {
        // check left macroblock pair first
          if (mb_is_available(img->current_mb_nr-2, img->current_mb_nr)&&((img->current_mb_nr%(img->PicWidthInMbs*2))!=0))
          {
            currMB->mb_field = img->mb_data[img->current_mb_nr-2].mb_field;
          }
          else
          {
            // check top macroblock pair
            if (mb_is_available(img->current_mb_nr-2*img->PicWidthInMbs, img->current_mb_nr))
            {
              currMB->mb_field = img->mb_data[img->current_mb_nr-2*img->PicWidthInMbs].mb_field;
            }
            else
              currMB->mb_field = 0;
          }
        }
      }
    }
  }

  dec_picture->mb_field[img->current_mb_nr] = currMB->mb_field;

  img->siblock[img->mb_x][img->mb_y]=0;

/**根据宏块类型的不同采用不同的interpret_mb_mode()函数来对宏块进行解析，如I、P、B、sI宏块。然后判断该宏块是否为
8x8子宏块分割，如果是，则解析子宏块分割模式。*********/
  if ((img->type==P_SLICE ))    // inter frame
    interpret_mb_mode_P(img);   //解码器出于程序的考虑将解码得到的宏块类型要进行一个转换，就是在这些 interpret 函数里实现的
  else if (img->type==I_SLICE)		// intra frame
    interpret_mb_mode_I(img);
  else if ((img->type==B_SLICE))       // B frame
    interpret_mb_mode_B(img);
  else if ((img->type==SP_SLICE))     // SP frame
    interpret_mb_mode_P(img);
  else if (img->type==SI_SLICE)     // SI frame
    interpret_mb_mode_SI(img);

  if(img->MbaffFrameFlag)
  {
    if(currMB->mb_field)
    {
      img->num_ref_idx_l0_active <<=1;
      img->num_ref_idx_l1_active <<=1;
    }
  }

  //====== READ 8x8 SUB-PARTITION MODES (modes of 8x8 blocks) and Intra VBST block modes ======
  //-------判断该宏块是否为8x8子宏块分割，如果是，则解析子宏块分割模式--------------------------
  if (IS_P8x8 (currMB))
  {
    currSE.type    = SE_MBTYPE;
    dP = &(currSlice->partArr[partMap[SE_MBTYPE]]);

    for (i=0; i<4; i++)
    {
      if (active_pps->entropy_coding_mode_flag ==UVLC || dP->bitstream->ei_flag) currSE.mapping = linfo_ue;
      else                                                  currSE.reading = readB8_typeInfo_CABAC;

      TRACE_STRING("sub_mb_type");
      dP->readSyntaxElement (&currSE, img, inp, dP);
      SetB8Mode (img, currMB, currSE.value1, i);	//++ 解析子宏块的分割模式
    }
  }

  if(active_pps->constrained_intra_pred_flag && (img->type==P_SLICE|| img->type==B_SLICE))        // inter frame
  {
    if( !IS_INTRA(currMB) )
    {
      img->intra_block[img->current_mb_nr] = 0;
    }
  }

  //! TO for Error Concelament
  //! If we have an INTRA Macroblock and we lost the partition
  //! which contains the intra coefficients Copy MB would be better 
  //! than just a grey block.
  //! Seems to be a bit at the wrong place to do this right here, but for this case 
  //! up to now there is no other way.
  //++ 为了进行错误隐藏，如果某个帧内宏块包含帧内系数的部分在传输过程中丢失了。则采用复制宏块的方法，而不是简单地用一个
  //++ 灰度值128的灰度块来填充。虽然这样看起来仍然不正确，但是对于这类宏块分割丢失情况，目前还没有其他处理方式。
    dP = &(currSlice->partArr[partMap[SE_CBP_INTRA]]);
  if(IS_INTRA (currMB) && dP->bitstream->ei_flag && img->number)	//++ 如果当前宏块是帧内块而且有错误而且不是第一个宏块，则
  {																	//++ 对当前宏块做一些设置，以便后面采用上面提到的复制宏块
    currMB->mb_type = 0;											//++ 的处理方法。另：因为第一个宏块没有可复制的对象，所以
    currMB->ei_flag = 1;											//++ 这里要判断当前宏块是不是第一个宏块。
    for (i=0;i<4;i++) {currMB->b8mode[i]=currMB->b8pdir[i]=0; }
  }
  dP = &(currSlice->partArr[partMap[currSE.type]]);
  //! End TO


  // init macroblock data ---init_macroblock()函数主要作用是对亮度块中16个4x4块的预测模式、运动矢量进行初始化。----
  init_macroblock       (img);	//++ 将亮度块中的16个4*4块的预测模式设置为2（直接预测），运动向量置0

  //---IS_DIRECT判断该宏块是否B_kip类型宏块，如果是，则按相应算法处理。接着IS_COPY判断该宏块是否COPY类型宏块，
  //如果是，则按相应算法处理。下一步，read_iPredes_modes从比特流中读取亮度块和色度块的预测模式。------------
  if (IS_DIRECT (currMB) && img->cod_counter >= 0)	//++ B_Skip类型宏块。无残差，无MVD。解码时，通过Direct预测模式（时间或空间）计算出前、后向MV（在 decode_one_macroblock 函数中计算），
  {													//++ 然后直接利用前、后向MV得到像素预测值。像素重构值=像素预测值
    currMB->cbp = 0;								//++ 如果img->cod_counter<0，则为B_Direct_16*16类型宏块。有残差，无MVD。解码过程与B_Skip类型宏块相同，只是多了一个残差的解码过程
    reset_coeffs();									//++ 即，像素重构值=像素预测值+残差解码值

    if (active_pps->entropy_coding_mode_flag ==CABAC)
      img->cod_counter=-1;

    return DECODE_MB;
  }

  if (IS_COPY (currMB)) //keep last macroblock	//++ COPY宏块（也即P_Skip类型宏块）。无残差，无MVD。直接利用预测MV得到像素预测值。像素重构值=像素预测值
  {
    int i, j, k, pmv[2];
    int zeroMotionAbove;
    int zeroMotionLeft;
    PixelPos mb_a, mb_b;
    int      a_mv_y = 0;
    int      a_ref_idx = 0;
    int      b_mv_y = 0;
    int      b_ref_idx = 0;
    int      list_offset = ((img->MbaffFrameFlag)&&(currMB->mb_field))? img->current_mb_nr%2 ? 4 : 2 : 0;

    getLuma4x4Neighbour(img->current_mb_nr,0,0,-1, 0,&mb_a);
    getLuma4x4Neighbour(img->current_mb_nr,0,0, 0,-1,&mb_b);

    if (mb_a.available)
    {
      a_mv_y    = dec_picture->mv[LIST_0][mb_a.pos_x][mb_a.pos_y][1];
      a_ref_idx = dec_picture->ref_idx[LIST_0][mb_a.pos_x][mb_a.pos_y];

      if (currMB->mb_field && !img->mb_data[mb_a.mb_addr].mb_field)
      {
        a_mv_y    /=2;
        a_ref_idx *=2;
      }
      if (!currMB->mb_field && img->mb_data[mb_a.mb_addr].mb_field)
      {
        a_mv_y    *=2;
        a_ref_idx >>=1;
      }
    }

    if (mb_b.available)
    {
      b_mv_y    = dec_picture->mv[LIST_0][mb_b.pos_x][mb_b.pos_y][1];
      b_ref_idx = dec_picture->ref_idx[LIST_0][mb_b.pos_x][mb_b.pos_y];

      if (currMB->mb_field && !img->mb_data[mb_b.mb_addr].mb_field)
      {
        b_mv_y    /=2;
        b_ref_idx *=2;
      }
      if (!currMB->mb_field && img->mb_data[mb_b.mb_addr].mb_field)
      {
        b_mv_y    *=2;
        b_ref_idx >>=1;
      }
    }
    
    zeroMotionLeft  = !mb_a.available ? 1 : a_ref_idx==0 && dec_picture->mv[LIST_0][mb_a.pos_x][mb_a.pos_y][0]==0 && a_mv_y==0 ? 1 : 0;
    zeroMotionAbove = !mb_b.available ? 1 : b_ref_idx==0 && dec_picture->mv[LIST_0][mb_b.pos_x][mb_b.pos_y][0]==0 && b_mv_y==0 ? 1 : 0;

    currMB->cbp = 0;
    reset_coeffs();

    img_block_y   = img->block_y;
  
    if (zeroMotionAbove || zeroMotionLeft)
    {
      for(i=0;i<BLOCK_SIZE;i++)
        for(j=0;j<BLOCK_SIZE;j++)
          for (k=0;k<2;k++)
            dec_picture->mv[LIST_0][img->block_x+i][img->block_y+j][k] = 0;
    }
    else
    {
      SetMotionVectorPredictor (img, pmv, pmv+1, 0, LIST_0, dec_picture->ref_idx, dec_picture->mv, 0, 0, 16, 16);
      
      for(i=0;i<BLOCK_SIZE;i++)
        for(j=0;j<BLOCK_SIZE;j++)
          for (k=0;k<2;k++)
          {
            dec_picture->mv[LIST_0][img->block_x+i][img_block_y+j][k] = pmv[k];
          }
    }

    for(i=0;i<BLOCK_SIZE;i++)
      for(j=0;j<BLOCK_SIZE;j++)
      {
        dec_picture->ref_idx[LIST_0][img->block_x+i][img_block_y+j] = 0;
        dec_picture->ref_pic_id[LIST_0][img->block_x+i][img_block_y+j] = dec_picture->ref_pic_num[img->current_slice_nr][LIST_0 + list_offset][dec_picture->ref_idx[LIST_0][img->block_x+i][img_block_y+j]];
      }

    return DECODE_MB;
  }
  if(currMB->mb_type!=IPCM)
  {
    
    // intra prediction modes for a macroblock 4x4 **********************************************
    read_ipred_modes(img,inp);	//++ 从比特流中读取亮度块和色度块的帧内预测模式,以亮度块为例，读取时是以4x4块为单位
	                            //读取的（有一个读取顺序，见标准）。
    // read inter frame vector data *********************************************************
    if (IS_INTERMV (currMB)) //---IS_INTERMV 判断该宏块是否帧间编码宏块，如果是，则通过函数readMotionInfoFromNAl()
		                     //从码流中读取运动矢量残差。---如果当前宏块不是帧内宏块、I_PCM宏块、B_Direct_16*16和
    {						//++ P_Skip、B_Skip宏块则从码流中读取运动向量差值
      readMotionInfoFromNAL (img, inp);//该函数分别进行了参考帧、MV的读取，这正是进行了 7.3.5.1 所规定的内容 
    }
    // read CBP and Coeffs  ***************************************************************
    readCBPandCoeffsFromNAL (img,inp);	//++ 从比特流中读取CBP和变换系数（包括反量化）
  }
  else
  {
    //read pcm_alignment_zero_bit and pcm_byte[i] 
    
    // here dP is assigned with the same dP as SE_MBTYPE, because IPCM syntax is in the 
    // same category as MBTYPE
    dP = &(currSlice->partArr[partMap[SE_MBTYPE]]);
    readIPCMcoeffsFromNAL(img,inp,dP);
  }
  return DECODE_MB;
}



/*!
 ************************************************************************
 * \brief
 *    Initialize decoding engine after decoding an IPCM macroblock
 *    (for IPCM CABAC  28/11/2003)
 *
 * \author
 *    Dong Wang <Dong.Wang@bristol.ac.uk>  
 ************************************************************************
 */
void init_decoding_engine_IPCM(struct img_par *img)
{
  Slice *currSlice = img->currentSlice;
  Bitstream *currStream;
  int ByteStartPosition;
  int PartitionNumber;
  int i;
  
  if(currSlice->dp_mode==PAR_DP_1)
    PartitionNumber=1;
  else if(currSlice->dp_mode==PAR_DP_3)
    PartitionNumber=3;
  else
  {
    printf("Partition Mode is not supported\n");
    exit(1);
  }
  
  for(i=0;i<PartitionNumber;i++)
  {
    currStream = currSlice->partArr[i].bitstream;
    ByteStartPosition = currStream->read_len;
    
    
    arideco_start_decoding (&currSlice->partArr[i].de_cabac, currStream->streamBuffer, ByteStartPosition, &currStream->read_len, img->type);
  }
}




/*!
 ************************************************************************
 * \brief
 *    Read IPCM pcm_alignment_zero_bit and pcm_byte[i] from stream to img->cof
 *    (for IPCM CABAC and IPCM CAVLC  28/11/2003)
 *
 * \author
 *    Dong Wang <Dong.Wang@bristol.ac.uk>
 ************************************************************************
 */

void readIPCMcoeffsFromNAL(struct img_par *img, struct inp_par *inp, struct datapartition *dP)
{
  SyntaxElement currSE;
  int i,j;
  
  //For CABAC, we don't need to read bits to let stream byte aligned
  //  because we have variable for integer bytes position
  if(active_pps->entropy_coding_mode_flag  == CABAC)
  {
    //read luma and chroma IPCM coefficients
    currSE.len=8;

  for(i=0;i<16;i++)
    for(j=0;j<16;j++)
    {
      readIPCMBytes_CABAC(&currSE, dP->bitstream);
      img->cof[i/4][j/4][i%4][j%4]=currSE.value1;
    }

  for(i=0;i<8;i++)
    for(j=0;j<8;j++)
    {
      readIPCMBytes_CABAC(&currSE, dP->bitstream);
      img->cof[i/4][j/4+4][i%4][j%4]=currSE.value1;
    }

  for(i=0;i<8;i++)
    for(j=0;j<8;j++)
    {
      readIPCMBytes_CABAC(&currSE, dP->bitstream);
      img->cof[i/4+2][j/4+4][i%4][j%4]=currSE.value1;
    }

  //If the decoded MB is IPCM MB, decoding engine is initialized

  // here the decoding engine is directly initialized without checking End of Slice
  // The reason is that, whether current MB is the last MB in slice or not, there is
  // at least one 'end of slice' syntax after this MB. So when fetching bytes in this 
    // initialisation process, we can guarantee there is bits available in bitstream. 
  init_decoding_engine_IPCM(img);

  }
  else
  { 

  //read bits to let stream byte aligned


  if((dP->bitstream->frame_bitoffset)%8!=0)
  {

    currSE.len=8-(dP->bitstream->frame_bitoffset)%8;
    readSyntaxElement_FLC(&currSE, dP->bitstream);
  
  }

  //read luma and chroma IPCM coefficients
  currSE.len=8;

  for(i=0;i<16;i++)
    for(j=0;j<16;j++)
    {
      readSyntaxElement_FLC(&currSE, dP->bitstream);
      img->cof[i/4][j/4][i%4][j%4]=currSE.value1;
    }

  for(i=0;i<8;i++)
    for(j=0;j<8;j++)
    {
      readSyntaxElement_FLC(&currSE, dP->bitstream);
      img->cof[i/4][j/4+4][i%4][j%4]=currSE.value1;
    }

  for(i=0;i<8;i++)
    for(j=0;j<8;j++)
    {
      readSyntaxElement_FLC(&currSE, dP->bitstream);
      img->cof[i/4+2][j/4+4][i%4][j%4]=currSE.value1;
    }
  }
}


//++ 从比特流中读取亮度块和色度块的预测模式。亮度块读取顺序如下：为了解谜做了一些结构变动
//++ ***************
//++ * 1  2  5  7  *
//++ * 3  4  6  8  *
//++ * 9  10 13 14 *
//++ * 11 12 15 16 *
//++ ***************
void read_ipred_modes(struct img_par *img,struct inp_par *inp)
{
  int b8,i,j,bi,bj,bx,by,dec;
  SyntaxElement currSE;
  Slice *currSlice;
  DataPartition *dP;
  int *partMap;
  Macroblock *currMB;
  int ts, ls;
  int mostProbableIntraPredMode;
  int upIntraPredMode;
  int leftIntraPredMode;
  int IntraChromaPredModeFlag;

  PixelPos left_block;
  PixelPos top_block;
  //////////////////////////////////////////////////////////////////
  int temp_intr;
  //////////////////////////////////////////////////////////////////
  currMB = &img->mb_data[img->current_mb_nr];

  IntraChromaPredModeFlag = IS_INTRA(currMB);

  currSlice = img->currentSlice;
  partMap = assignSE2partition[currSlice->dp_mode];

  currSE.type = SE_INTRAPREDMODE;

  TRACE_STRING("intra4x4_pred_mode");
  dP = &(currSlice->partArr[partMap[currSE.type]]);

  if (!(active_pps->entropy_coding_mode_flag == UVLC || dP->bitstream->ei_flag)) 
    currSE.reading = readIntraPredMode_CABAC;

  for(b8=0;b8<4;b8++)  //loop 8x8 blocks	//++ 该循环读取4个8*8的块
  {
    if( currMB->b8mode[b8]==IBLOCK )
    {
      IntraChromaPredModeFlag = 1;

      for(j=0;j<2;j++)  //loop subblocks	//++ 该双重循环读取亮度块中8*8块中4个4*4块的预测模式
        for(i=0;i<2;i++)
        {
          //get from stream
          if (active_pps->entropy_coding_mode_flag == UVLC || dP->bitstream->ei_flag)
            readSyntaxElement_Intra4x4PredictionMode(&currSE,img,inp,dP);
          else 
          {
            currSE.context=(b8<<2)+(j<<1)+i;
            dP->readSyntaxElement(&currSE,img,inp,dP);
          }

          bx = ((b8&1)<<1) + i;
          by = (b8&2)      + j;

          //getLuma4x4Neighbour(img->current_mb_nr, bx, by, -1,  0, &left_block);
          //getLuma4x4Neighbour(img->current_mb_nr, bx, by,  0, -1, &top_block);
          
          //get from array and decode
          bi = img->block_x + bx;
          bj = img->block_y + by;

         /* if (active_pps->constrained_intra_pred_flag)
          {
            left_block.available = left_block.available ? img->intra_block[left_block.mb_addr] : 0;
            top_block.available  = top_block.available  ? img->intra_block[top_block.mb_addr]  : 0;
          }

          // !! KS: not sure if the follwing is still correct...
          ts=ls=0;   // Check to see if the neighboring block is SI
          if (IS_OLDINTRA(currMB) && img->type == SI_SLICE)           // need support for MBINTLC1
          {
            if (left_block.available)
              if (img->siblock [left_block.pos_x][left_block.pos_y])
                ls=1;

            if (top_block.available)
              if (img->siblock [top_block.pos_x][top_block.pos_y])
                ts=1;
          }

          upIntraPredMode            = (top_block.available  &&(ts == 0)) ? img->ipredmode[top_block.pos_x ][top_block.pos_y ] : -1;
          leftIntraPredMode          = (left_block.available &&(ls == 0)) ? img->ipredmode[left_block.pos_x][left_block.pos_y] : -1;

          mostProbableIntraPredMode  = (upIntraPredMode < 0 || leftIntraPredMode < 0) ? DC_PRED : upIntraPredMode < leftIntraPredMode ? upIntraPredMode : leftIntraPredMode;*/
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
          //dec = (currSE.value1 == -1) ? mostProbableIntraPredMode : currSE.value1 + (currSE.value1 >= mostProbableIntraPredMode);
          dec = currSE.value1; 
          //set
          img->ipredmode[bi][bj]=dec;
          }
        }
    }
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
if(en_ps.intra_mode)
{
 if(currMB->mb_type==I4MB)
  decryption_intra_mode(img);
}
for(b8=0;b8<4;b8++)  //loop 8x8 blocks	//++ 该循环读取4个8*8的块
  {
    if( currMB->b8mode[b8]==IBLOCK )
    {
      for(j=0;j<2;j++)  //loop subblocks	//++ 该双重循环读取亮度块中8*8块中4个4*4块的预测模式
        for(i=0;i<2;i++)
        {
         bx = ((b8&1)<<1) + i;
         by = (b8&2) + j;
         getLuma4x4Neighbour(img->current_mb_nr, bx, by, -1,  0, &left_block);
         getLuma4x4Neighbour(img->current_mb_nr, bx, by,  0, -1, &top_block);
         bi = img->block_x + bx;
         bj = img->block_y + by;
         if (active_pps->constrained_intra_pred_flag)
         {
            left_block.available = left_block.available ? img->intra_block[left_block.mb_addr] : 0;
            top_block.available  = top_block.available  ? img->intra_block[top_block.mb_addr]  : 0;
         }
         // !! KS: not sure if the follwing is still correct...
         ts=ls=0;// Check to see if the neighboring block is SI
         if (IS_OLDINTRA(currMB) && img->type == SI_SLICE)           // need support for MBINTLC1
          {
            if (left_block.available)
              if (img->siblock [left_block.pos_x][left_block.pos_y])
                ls=1;

            if (top_block.available)
              if (img->siblock [top_block.pos_x][top_block.pos_y])
                ts=1;
          }
          upIntraPredMode            = (top_block.available  &&(ts == 0)) ? img->ipredmode[top_block.pos_x ][top_block.pos_y ] : -1;
          leftIntraPredMode          = (left_block.available &&(ls == 0)) ? img->ipredmode[left_block.pos_x][left_block.pos_y] : -1;
          mostProbableIntraPredMode  = (upIntraPredMode < 0 || leftIntraPredMode < 0) ? DC_PRED : upIntraPredMode < leftIntraPredMode ? upIntraPredMode : leftIntraPredMode;
          temp_intr= img->ipredmode[bi][bj];
          dec=(temp_intr==-1) ? mostProbableIntraPredMode : temp_intr+(temp_intr>=mostProbableIntraPredMode);
          img->ipredmode[bi][bj]=dec;
         }
     }
 }
  if (IntraChromaPredModeFlag)
  {
    currSE.type = SE_INTRAPREDMODE;
    TRACE_STRING("intra_chroma_pred_mode");
    dP = &(currSlice->partArr[partMap[currSE.type]]);

    if (active_pps->entropy_coding_mode_flag == UVLC || dP->bitstream->ei_flag) currSE.mapping = linfo_ue;
    else                                                    currSE.reading = readCIPredMode_CABAC;

    dP->readSyntaxElement(&currSE,img,inp,dP);
    currMB->c_ipred_mode = currSE.value1;
    ///////////////////////////////////////////////////////////////////////////////////////////////
    if(en_ps.chroma_mode)
     currMB->c_ipred_mode = decryption_C_intra_mode(currMB->c_ipred_mode);
    ///////////////////////////////////////////////////////////////////////////////////////////////
    if (currMB->c_ipred_mode < DC_PRED_8 || currMB->c_ipred_mode > PLANE_8)
    {
      error("illegal chroma intra pred mode!\n", 600);
    }
  }
}



/*!
 ************************************************************************
 * \brief
 *    Set motion vector predictor
 ************************************************************************
 */
static void SetMotionVectorPredictor (struct img_par  *img,
                                      int             *pmv_x,
                                      int             *pmv_y,
                                      int             ref_frame,
                                      int             list,
                                      int             ***refPic,
                                      int             ****tmp_mv,
                                      int             block_x,	//++ 当前分割左上角 4*4 块在本宏块内的横坐标（以 4*4 块为单位）
                                      int             block_y,	//++ 当前分割左上角 4*4 块在本宏块内的纵坐标（以 4*4 块为单位）
                                      int             blockshape_x,
                                      int             blockshape_y)
{	//++ 该函数计算运动向量预测值
  int mb_x                 = BLOCK_SIZE*block_x;	//++ 当前分割左上角像素点在本宏块内的横坐标
  int mb_y                 = BLOCK_SIZE*block_y;	//++ 当前分割左上角像素点在本宏块内的纵坐标
  int mb_nr                = img->current_mb_nr;

  int mv_a, mv_b, mv_c, pred_vec=0;
  int mvPredType, rFrameL, rFrameU, rFrameUR;
  int hv;


  PixelPos block_a, block_b, block_c, block_d;

  //++ 下面的四个 getLuma4x4Neighbour 函数只是判断 A、B、C、D 所在宏块是否可用
  getLuma4x4Neighbour(mb_nr, block_x, block_y,           -1,  0, &block_a);
  getLuma4x4Neighbour(mb_nr, block_x, block_y,            0, -1, &block_b);
  getLuma4x4Neighbour(mb_nr, block_x, block_y, blockshape_x, -1, &block_c);
  getLuma4x4Neighbour(mb_nr, block_x, block_y,           -1, -1, &block_d);

  //++ 当 A、B、C、D 所在宏块为本宏块时（即 A、B、C、D 在本宏块内部），则需要进一步判断他们的可用性（因为在内部并不一定可用）
  if (mb_y > 0)
  {
    if (mb_x < 8)  // first column of 8x8 blocks
    {
      if (mb_y==8)
      {
        if (blockshape_x == 16)      block_c.available  = 0;
        else                         block_c.available &= 1;
      }
      else
      {
        if (mb_x+blockshape_x != 8)  block_c.available &= 1;
        else                         block_c.available  = 0;
      }
    }
    else
    {
      if (mb_x+blockshape_x != 16)   block_c.available &= 1;
      else                           block_c.available  = 0;
    }
  }

  if (!block_c.available)
  {
    block_c=block_d;
  }

  mvPredType = MVPRED_MEDIAN;

  if (!img->MbaffFrameFlag)
  {
    rFrameL    = block_a.available    ? refPic[list][block_a.pos_x][block_a.pos_y] : -1;
    rFrameU    = block_b.available    ? refPic[list][block_b.pos_x][block_b.pos_y] : -1;
    rFrameUR   = block_c.available    ? refPic[list][block_c.pos_x][block_c.pos_y] : -1;
  }
  else
  {
    if (img->mb_data[img->current_mb_nr].mb_field)
    {
      rFrameL    = block_a.available    ? 
        img->mb_data[block_a.mb_addr].mb_field ? 
        refPic[list][block_a.pos_x][block_a.pos_y]:
        refPic[list][block_a.pos_x][block_a.pos_y] * 2: 
        -1;
      rFrameU    = block_b.available    ? 
        img->mb_data[block_b.mb_addr].mb_field ? 
        refPic[list][block_b.pos_x][block_b.pos_y]:
        refPic[list][block_b.pos_x][block_b.pos_y] * 2: 
        -1;
      rFrameUR    = block_c.available    ? 
        img->mb_data[block_c.mb_addr].mb_field ? 
        refPic[list][block_c.pos_x][block_c.pos_y]:
        refPic[list][block_c.pos_x][block_c.pos_y] * 2: 
        -1;
    }
    else
    {
      rFrameL    = block_a.available    ? 
        img->mb_data[block_a.mb_addr].mb_field ? 
        refPic[list][block_a.pos_x][block_a.pos_y] >>1:
        refPic[list][block_a.pos_x][block_a.pos_y] : 
        -1;
      rFrameU    = block_b.available    ? 
        img->mb_data[block_b.mb_addr].mb_field ? 
        refPic[list][block_b.pos_x][block_b.pos_y] >>1:
        refPic[list][block_b.pos_x][block_b.pos_y] : 
        -1;
      rFrameUR    = block_c.available    ? 
        img->mb_data[block_c.mb_addr].mb_field ? 
        refPic[list][block_c.pos_x][block_c.pos_y] >>1:
        refPic[list][block_c.pos_x][block_c.pos_y] : 
        -1;
    }
  }


  /* Prediction if only one of the neighbors uses the reference frame
   * we are checking
   */
  if(rFrameL == ref_frame && rFrameU != ref_frame && rFrameUR != ref_frame)       mvPredType = MVPRED_L;
  else if(rFrameL != ref_frame && rFrameU == ref_frame && rFrameUR != ref_frame)  mvPredType = MVPRED_U;
  else if(rFrameL != ref_frame && rFrameU != ref_frame && rFrameUR == ref_frame)  mvPredType = MVPRED_UR;
  // Directional predictions 
  if(blockshape_x == 8 && blockshape_y == 16)
  {
    if(mb_x == 0)
    {
      if(rFrameL == ref_frame)
        mvPredType = MVPRED_L;
    }
    else
    {
      if( rFrameUR == ref_frame)
        mvPredType = MVPRED_UR;
    }
  }
  else if(blockshape_x == 16 && blockshape_y == 8)
  {
    if(mb_y == 0)
    {
      if(rFrameU == ref_frame)
        mvPredType = MVPRED_U;
    }
    else
    {
      if(rFrameL == ref_frame)
        mvPredType = MVPRED_L;
    }
  }

  for (hv=0; hv < 2; hv++)
  {
    if (!img->MbaffFrameFlag || hv==0)
    {
      mv_a = block_a.available  ? tmp_mv[list][block_a.pos_x][block_a.pos_y][hv] : 0;
      mv_b = block_b.available  ? tmp_mv[list][block_b.pos_x][block_b.pos_y][hv] : 0;
      mv_c = block_c.available  ? tmp_mv[list][block_c.pos_x][block_c.pos_y][hv] : 0;
    }
    else
    {
      if (img->mb_data[img->current_mb_nr].mb_field)
      {
        mv_a = block_a.available  ? img->mb_data[block_a.mb_addr].mb_field?
          tmp_mv[list][block_a.pos_x][block_a.pos_y][hv]:
          tmp_mv[list][block_a.pos_x][block_a.pos_y][hv] / 2: 
          0;
        mv_b = block_b.available  ? img->mb_data[block_b.mb_addr].mb_field?
          tmp_mv[list][block_b.pos_x][block_b.pos_y][hv]:
          tmp_mv[list][block_b.pos_x][block_b.pos_y][hv] / 2: 
          0;
        mv_c = block_c.available  ? img->mb_data[block_c.mb_addr].mb_field?
          tmp_mv[list][block_c.pos_x][block_c.pos_y][hv]:
          tmp_mv[list][block_c.pos_x][block_c.pos_y][hv] / 2: 
          0;
      }
      else
      {
        mv_a = block_a.available  ? img->mb_data[block_a.mb_addr].mb_field?
          tmp_mv[list][block_a.pos_x][block_a.pos_y][hv] * 2:
          tmp_mv[list][block_a.pos_x][block_a.pos_y][hv]: 
          0;
        mv_b = block_b.available  ? img->mb_data[block_b.mb_addr].mb_field?
          tmp_mv[list][block_b.pos_x][block_b.pos_y][hv] * 2:
          tmp_mv[list][block_b.pos_x][block_b.pos_y][hv]: 
          0;
        mv_c = block_c.available  ? img->mb_data[block_c.mb_addr].mb_field?
          tmp_mv[list][block_c.pos_x][block_c.pos_y][hv] * 2:
          tmp_mv[list][block_c.pos_x][block_c.pos_y][hv]: 
          0;
      }
    }


    switch (mvPredType)
    {
    case MVPRED_MEDIAN:
      if(!(block_b.available || block_c.available))
        pred_vec = mv_a;
      else
        pred_vec = mv_a+mv_b+mv_c-min(mv_a,min(mv_b,mv_c))-max(mv_a,max(mv_b,mv_c));
      break;
    case MVPRED_L:
      pred_vec = mv_a;
      break;
    case MVPRED_U:
      pred_vec = mv_b;
      break;
    case MVPRED_UR:
      pred_vec = mv_c;
      break;
    default:
      break;
    }

    if (hv==0)  *pmv_x = pred_vec;
    else        *pmv_y = pred_vec;

  }
}


/*!
 ************************************************************************
 * \brief
 *    Set context for reference frames
 ************************************************************************
 */
int
BType2CtxRef (int btype)
{
  if (btype<4)  return 0;
  else          return 1;
}

/*!
 ************************************************************************
 * \brief
 *    Read motion info
 ************************************************************************
 */
void readMotionInfoFromNAL (struct img_par *img, struct inp_par *inp)
{
  int i,j,k;
  int step_h,step_v;
  int curr_mvd;
  Macroblock *currMB  = &img->mb_data[img->current_mb_nr];
  SyntaxElement currSE;
  Slice *currSlice    = img->currentSlice;
  DataPartition *dP;
  int *partMap        = assignSE2partition[currSlice->dp_mode];
  int bframe          = (img->type==B_SLICE);
  int partmode        = (IS_P8x8(currMB)?4:currMB->mb_type);
  int step_h0         = BLOCK_STEP [partmode][0];
  int step_v0         = BLOCK_STEP [partmode][1];

  int mv_mode, i0, j0, refframe;
  int pmv[2];
  int j4, i4, ii,jj;
  int vec;

  int mv_scale = 0;

  int flag_mode;

  int list_offset = ((img->MbaffFrameFlag)&&(currMB->mb_field))? img->current_mb_nr%2 ? 4 : 2 : 0;

  byte **    moving_block;
  int ****   co_located_mv;
  int ***    co_located_ref_idx;
  int64 ***    co_located_ref_id;

  if ((img->MbaffFrameFlag)&&(currMB->mb_field))
  {
    if(img->current_mb_nr%2)
    {
      moving_block = Co_located->bottom_moving_block;
      co_located_mv = Co_located->bottom_mv;
      co_located_ref_idx = Co_located->bottom_ref_idx;
      co_located_ref_id = Co_located->bottom_ref_pic_id;
    }
    else
    {
      moving_block = Co_located->top_moving_block;
      co_located_mv = Co_located->top_mv;
      co_located_ref_idx = Co_located->top_ref_idx;
      co_located_ref_id = Co_located->top_ref_pic_id;
    }
  }
  else
  {
    moving_block = Co_located->moving_block;
    co_located_mv = Co_located->mv;
    co_located_ref_idx = Co_located->ref_idx;
    co_located_ref_id = Co_located->ref_pic_id;
  }

  if (bframe && IS_P8x8 (currMB))	//++ 计算B_Direct_8*8类型子块的参考帧
  {
    if (img->direct_type)	//++ 空间直接模式，这里主要是因为可能子宏块分割中存在Direct模式，所以先进行一个时空Direct模式的判断
    {
      int imgblock_y= ((img->MbaffFrameFlag)&&(currMB->mb_field))? (img->current_mb_nr%2) ? (img->block_y-4)/2:img->block_y/2: img->block_y;
      int fw_rFrameL, fw_rFrameU, fw_rFrameUL, fw_rFrameUR;
      int bw_rFrameL, bw_rFrameU, bw_rFrameUL, bw_rFrameUR;
      
      PixelPos mb_left, mb_up, mb_upleft, mb_upright;
      
      int fw_rFrame,bw_rFrame;
      int pmvfw[2]={0,0},pmvbw[2]={0,0};
     
      
      getLuma4x4Neighbour(img->current_mb_nr, 0, 0, -1,  0, &mb_left);
      getLuma4x4Neighbour(img->current_mb_nr, 0, 0,  0, -1, &mb_up);
      getLuma4x4Neighbour(img->current_mb_nr, 0, 0, 16, -1, &mb_upright);
      getLuma4x4Neighbour(img->current_mb_nr, 0, 0, -1, -1, &mb_upleft);

      if (!img->MbaffFrameFlag)
      {
        fw_rFrameL = mb_left.available ? dec_picture->ref_idx[LIST_0][mb_left.pos_x][mb_left.pos_y] : -1;
        fw_rFrameU = mb_up.available ? dec_picture->ref_idx[LIST_0][mb_up.pos_x][mb_up.pos_y] : -1;
        fw_rFrameUL = mb_upleft.available ? dec_picture->ref_idx[LIST_0][mb_upleft.pos_x][mb_upleft.pos_y] : -1;
        fw_rFrameUR = mb_upright.available ? dec_picture->ref_idx[LIST_0][mb_upright.pos_x][mb_upright.pos_y] : fw_rFrameUL;      
        
        bw_rFrameL = mb_left.available ? dec_picture->ref_idx[LIST_1][mb_left.pos_x][mb_left.pos_y] : -1;
        bw_rFrameU = mb_up.available ? dec_picture->ref_idx[LIST_1][mb_up.pos_x][mb_up.pos_y] : -1;
        bw_rFrameUL = mb_upleft.available ? dec_picture->ref_idx[LIST_1][mb_upleft.pos_x][mb_upleft.pos_y] : -1;
        bw_rFrameUR = mb_upright.available ? dec_picture->ref_idx[LIST_1][mb_upright.pos_x][mb_upright.pos_y] : bw_rFrameUL;      
      }
      else
      {
        if (img->mb_data[img->current_mb_nr].mb_field)
        {
          fw_rFrameL = mb_left.available ? 
            img->mb_data[mb_left.mb_addr].mb_field  || dec_picture->ref_idx[LIST_0][mb_left.pos_x][mb_left.pos_y] < 0? 
            dec_picture->ref_idx[LIST_0][mb_left.pos_x][mb_left.pos_y] : 
          dec_picture->ref_idx[LIST_0][mb_left.pos_x][mb_left.pos_y] * 2: -1;
          fw_rFrameU = mb_up.available ? 
            img->mb_data[mb_up.mb_addr].mb_field || dec_picture->ref_idx[LIST_0][mb_up.pos_x][mb_up.pos_y] < 0? 
            dec_picture->ref_idx[LIST_0][mb_up.pos_x][mb_up.pos_y] : 
          dec_picture->ref_idx[LIST_0][mb_up.pos_x][mb_up.pos_y] * 2: -1;

          fw_rFrameUL = mb_upleft.available ? 
            img->mb_data[mb_upleft.mb_addr].mb_field || dec_picture->ref_idx[LIST_0][mb_upleft.pos_x][mb_upleft.pos_y] < 0?
            dec_picture->ref_idx[LIST_0][mb_upleft.pos_x][mb_upleft.pos_y] : 
          dec_picture->ref_idx[LIST_0][mb_upleft.pos_x][mb_upleft.pos_y] *2: -1;      

          fw_rFrameUR = mb_upright.available ? 
            img->mb_data[mb_upright.mb_addr].mb_field || dec_picture->ref_idx[LIST_0][mb_upright.pos_x][mb_upright.pos_y] < 0 ?
            dec_picture->ref_idx[LIST_0][mb_upright.pos_x][mb_upright.pos_y] : 
          dec_picture->ref_idx[LIST_0][mb_upright.pos_x][mb_upright.pos_y] * 2: fw_rFrameUL;      
          
          bw_rFrameL = mb_left.available ? 
            img->mb_data[mb_left.mb_addr].mb_field || dec_picture->ref_idx[LIST_1][mb_left.pos_x][mb_left.pos_y]  < 0 ?
            dec_picture->ref_idx[LIST_1][mb_left.pos_x][mb_left.pos_y] : 
          dec_picture->ref_idx[LIST_1][mb_left.pos_x][mb_left.pos_y] * 2: -1;

          bw_rFrameU = mb_up.available ? 
            img->mb_data[mb_up.mb_addr].mb_field || dec_picture->ref_idx[LIST_1][mb_up.pos_x][mb_up.pos_y]  < 0 ?
            dec_picture->ref_idx[LIST_1][mb_up.pos_x][mb_up.pos_y] : 
          dec_picture->ref_idx[LIST_1][mb_up.pos_x][mb_up.pos_y] * 2: -1;

          bw_rFrameUL = mb_upleft.available ? 
            img->mb_data[mb_upleft.mb_addr].mb_field || dec_picture->ref_idx[LIST_1][mb_upleft.pos_x][mb_upleft.pos_y]  < 0 ?
            dec_picture->ref_idx[LIST_1][mb_upleft.pos_x][mb_upleft.pos_y] : 
          dec_picture->ref_idx[LIST_1][mb_upleft.pos_x][mb_upleft.pos_y] *2: -1;      

          bw_rFrameUR = mb_upright.available ? 
            img->mb_data[mb_upright.mb_addr].mb_field || dec_picture->ref_idx[LIST_1][mb_upright.pos_x][mb_upright.pos_y]  < 0 ?
            dec_picture->ref_idx[LIST_1][mb_upright.pos_x][mb_upright.pos_y] : 
          dec_picture->ref_idx[LIST_1][mb_upright.pos_x][mb_upright.pos_y] * 2: bw_rFrameUL;      
          
        }
        else
        {
          fw_rFrameL = mb_left.available ? 
            img->mb_data[mb_left.mb_addr].mb_field || dec_picture->ref_idx[LIST_0][mb_left.pos_x][mb_left.pos_y]  < 0 ?
            dec_picture->ref_idx[LIST_0][mb_left.pos_x][mb_left.pos_y] >> 1 : 
          dec_picture->ref_idx[LIST_0][mb_left.pos_x][mb_left.pos_y]: -1;

          fw_rFrameU = mb_up.available ? 
            img->mb_data[mb_up.mb_addr].mb_field || dec_picture->ref_idx[LIST_0][mb_up.pos_x][mb_up.pos_y]  < 0 ?
            dec_picture->ref_idx[LIST_0][mb_up.pos_x][mb_up.pos_y] >> 1 :  
          dec_picture->ref_idx[LIST_0][mb_up.pos_x][mb_up.pos_y] : -1;

          fw_rFrameUL = mb_upleft.available ? 
            img->mb_data[mb_upleft.mb_addr].mb_field || dec_picture->ref_idx[LIST_0][mb_upleft.pos_x][mb_upleft.pos_y] < 0 ?
            dec_picture->ref_idx[LIST_0][mb_upleft.pos_x][mb_upleft.pos_y]>> 1 : 
          dec_picture->ref_idx[LIST_0][mb_upleft.pos_x][mb_upleft.pos_y] : -1;      

          fw_rFrameUR = mb_upright.available ? 
            img->mb_data[mb_upright.mb_addr].mb_field || dec_picture->ref_idx[LIST_0][mb_upright.pos_x][mb_upright.pos_y]  < 0 ?
            dec_picture->ref_idx[LIST_0][mb_upright.pos_x][mb_upright.pos_y] >> 1 :  
          dec_picture->ref_idx[LIST_0][mb_upright.pos_x][mb_upright.pos_y] : fw_rFrameUL;      
          
          bw_rFrameL = mb_left.available ? 
            img->mb_data[mb_left.mb_addr].mb_field || dec_picture->ref_idx[LIST_1][mb_left.pos_x][mb_left.pos_y] < 0 ?
            dec_picture->ref_idx[LIST_1][mb_left.pos_x][mb_left.pos_y] >> 1 :  
          dec_picture->ref_idx[LIST_1][mb_left.pos_x][mb_left.pos_y] : -1;
          bw_rFrameU = mb_up.available ? 
            img->mb_data[mb_up.mb_addr].mb_field || dec_picture->ref_idx[LIST_1][mb_up.pos_x][mb_up.pos_y] < 0 ?
            dec_picture->ref_idx[LIST_1][mb_up.pos_x][mb_up.pos_y] >> 1 : 
          dec_picture->ref_idx[LIST_1][mb_up.pos_x][mb_up.pos_y] : -1;

          bw_rFrameUL = mb_upleft.available ? 
            img->mb_data[mb_upleft.mb_addr].mb_field || dec_picture->ref_idx[LIST_1][mb_upleft.pos_x][mb_upleft.pos_y]  < 0 ?
            dec_picture->ref_idx[LIST_1][mb_upleft.pos_x][mb_upleft.pos_y] >> 1 : 
          dec_picture->ref_idx[LIST_1][mb_upleft.pos_x][mb_upleft.pos_y] : -1;      

          bw_rFrameUR = mb_upright.available ? 
            img->mb_data[mb_upright.mb_addr].mb_field || dec_picture->ref_idx[LIST_1][mb_upright.pos_x][mb_upright.pos_y] < 0 ?
            dec_picture->ref_idx[LIST_1][mb_upright.pos_x][mb_upright.pos_y] >> 1: 
          dec_picture->ref_idx[LIST_1][mb_upright.pos_x][mb_upright.pos_y] : bw_rFrameUL;      
        }
      }
      
    fw_rFrame = (fw_rFrameL >= 0 && fw_rFrameU >= 0) ? min(fw_rFrameL,fw_rFrameU): max(fw_rFrameL,fw_rFrameU);
    fw_rFrame = (fw_rFrame >= 0 && fw_rFrameUR >= 0) ? min(fw_rFrame,fw_rFrameUR): max(fw_rFrame,fw_rFrameUR);
      
    bw_rFrame = (bw_rFrameL >= 0 && bw_rFrameU >= 0) ? min(bw_rFrameL,bw_rFrameU): max(bw_rFrameL,bw_rFrameU);
    bw_rFrame = (bw_rFrame >= 0 && bw_rFrameUR >= 0) ? min(bw_rFrame,bw_rFrameUR): max(bw_rFrame,bw_rFrameUR);
      
      
    if (fw_rFrame >=0)
        SetMotionVectorPredictor (img, pmvfw, pmvfw+1, fw_rFrame, LIST_0, dec_picture->ref_idx, dec_picture->mv, 0, 0, 16, 16);
      
    if (bw_rFrame >=0)
        SetMotionVectorPredictor (img, pmvbw, pmvbw+1, bw_rFrame, LIST_1, dec_picture->ref_idx, dec_picture->mv, 0, 0, 16, 16);
      
      
      for (i=0;i<4;i++)
      {
        if (currMB->b8mode[i] == 0)
          for(j=2*(i/2);j<2*(i/2)+2;j++)
            for(k=2*(i%2);k<2*(i%2)+2;k++)
            {
              int j6 = imgblock_y+j;
              j4 = img->block_y+j;
              i4 = img->block_x+k;
              
              
              if (fw_rFrame >= 0)
              {
                
                if  (!fw_rFrame  && ((!moving_block[i4][j6]) && (!listX[1+list_offset][0]->is_long_term)))
                {                    
                  dec_picture->mv  [LIST_0][i4][j4][0] = 0;
                  dec_picture->mv  [LIST_0][i4][j4][1] = 0;
                  dec_picture->ref_idx[LIST_0][i4][j4] = 0;                    
                }
                else
                {
                  
                  dec_picture->mv  [LIST_0][i4][j4][0] = pmvfw[0];
                  dec_picture->mv  [LIST_0][i4][j4][1] = pmvfw[1];
                  dec_picture->ref_idx[LIST_0][i4][j4] = fw_rFrame;
                }
              }
              else
              {
                dec_picture->mv  [LIST_0][i4][j4][0] = 0;
                dec_picture->mv  [LIST_0][i4][j4][1] = 0;
                dec_picture->ref_idx[LIST_0][i4][j4] = -1;
              }
              if (bw_rFrame >= 0)
              {
                if  (bw_rFrame==0 && ((!moving_block[i4][j6])&& (!listX[1+list_offset][0]->is_long_term)))
                {
                  dec_picture->mv  [LIST_1][i4][j4][0] = 0;
                  dec_picture->mv  [LIST_1][i4][j4][1] = 0;
                  dec_picture->ref_idx[LIST_1][i4][j4] = 0;
                }
                else
                {
                  dec_picture->mv  [LIST_1][i4][j4][0] = pmvbw[0];
                  dec_picture->mv  [LIST_1][i4][j4][1] = pmvbw[1];
                  dec_picture->ref_idx[LIST_1][i4][j4] = bw_rFrame;
                }
              }
              else
              {
                dec_picture->mv  [LIST_1][i4][j4][0] = 0;
                dec_picture->mv  [LIST_1][i4][j4][1] = 0;
                dec_picture->ref_idx[LIST_1][i4][j4] = -1;                               
              }
              
              if (fw_rFrame <0 && bw_rFrame <0)
              {
                dec_picture->ref_idx[LIST_0][i4][j4] = 0;
                dec_picture->ref_idx[LIST_1][i4][j4] = 0;                  
              }
            }
      }
    }
    else	//++ 以下代码是计算时间直接模式的参考帧，但是后面又做了一次计算，因此注释掉此处
    {/*
      for (i=0;i<4;i++)
      {
        if (currMB->b8mode[i] == 0)	//++ 8*8子宏块为B_Direct_8*8类型
        {

          for(j=2*(i/2);j<2*(i/2)+2;j++)
          {
            for(k=2*(i%2);k<2*(i%2)+2;k++)
            {
              
              int list_offset = ((img->MbaffFrameFlag)&&(currMB->mb_field))? img->current_mb_nr%2 ? 4 : 2 : 0;
              int imgblock_y= ((img->MbaffFrameFlag)&&(currMB->mb_field))? (img->current_mb_nr%2) ? (img->block_y-4)/2 : img->block_y/2 : img->block_y;
              int refList = co_located_ref_idx[LIST_0 ][img->block_x+k][imgblock_y+j]== -1 ? LIST_1 : LIST_0;	//++ 如果“相对位置块”的前向参考帧不可用则采用其后向参考帧
              int ref_idx = co_located_ref_idx[refList][img->block_x + k][imgblock_y + j];
              int mapped_idx=-1, iref;                             
 


              if (ref_idx == -1)
              {
                dec_picture->ref_idx [LIST_0][img->block_x + k][img->block_y + j] = 0;
                dec_picture->ref_idx [LIST_1][img->block_x + k][img->block_y + j] = 0;                
              }
              else
              {
                for (iref=0;iref<min(img->num_ref_idx_l0_active,listXsize[LIST_0 + list_offset]);iref++)
                {
                  if (dec_picture->ref_pic_num[img->current_slice_nr][LIST_0 + list_offset][iref]==co_located_ref_id[refList][img->block_x + k][imgblock_y + j])
                  {
                    mapped_idx=iref;
                    break;
                  }
                  else //! invalid index. Default to zero even though this case should not happen
                    mapped_idx=INVALIDINDEX;
                }
                if (INVALIDINDEX == mapped_idx)
                {
                  error("temporal direct error\ncolocated block has ref that is unavailable",-1111);
                }
                dec_picture->ref_idx [LIST_0][img->block_x + k][img->block_y + j] = mapped_idx;
                dec_picture->ref_idx [LIST_1][img->block_x + k][img->block_y + j] = 0;                
              }
            }
          }
        }
      }
  */}
  } 

  //  If multiple ref. frames, read reference frame for the MB *********************************
  if(img->num_ref_idx_l0_active>1)	//++ 计算非Direc类型块分割的前向参考帧（多帧参考）
  {
    flag_mode = ( img->num_ref_idx_l0_active == 2 ? 1 : 0);	//++ 如果参考帧数只有两帧，则采用一比特定长编码。如果参考帧数大于2，则采用变长编码（哥伦布编码）
															//++ 这是因为如果参考帧数为2，采用FLC只需要一个比特就能表示参考帧号，而VLC可能需要三个比特（当参考帧号不为0时）
    currSE.type = SE_REFFRAME;
    dP = &(currSlice->partArr[partMap[SE_REFFRAME]]);

    if (active_pps->entropy_coding_mode_flag == UVLC || dP->bitstream->ei_flag)   currSE.mapping = linfo_ue;
    else                                                      currSE.reading = readRefFrame_CABAC;
    
    for (j0=0; j0<4; j0+=step_v0)
    {
      for (i0=0; i0<4; i0+=step_h0)
      {
        k=2*(j0/2)+(i0/2);
        if ((currMB->b8pdir[k]==0 || currMB->b8pdir[k]==2) && currMB->b8mode[k]!=0)
        {
          TRACE_STRING("ref_idx_l0");

          img->subblock_x = i0;
          img->subblock_y = j0;
          
          if (!IS_P8x8 (currMB) || bframe || (!bframe && !img->allrefzero))
          {
            currSE.context = BType2CtxRef (currMB->b8mode[k]);	//++ 用于CABAC
            if( (active_pps->entropy_coding_mode_flag == UVLC || dP->bitstream->ei_flag) && flag_mode )
            {
              currSE.len = 1;
              readSyntaxElement_FLC(&currSE, dP->bitstream);
              currSE.value1 = 1 - currSE.value1;
            }
            else
            {
              currSE.value2 = LIST_0;
              dP->readSyntaxElement (&currSE,img,inp,dP);
            }
            refframe = currSE.value1;
            
          }
          else
          {
            refframe = 0;
          }
          
          /*
          if (bframe && refframe>img->buf_cycle)    // img->buf_cycle should be correct for field MBs now
          {
            set_ec_flag(SE_REFFRAME);
            refframe = 1;
          }
          */
          
          for (j=j0; j<j0+step_v0;j++)
            for (i=i0; i<i0+step_h0;i++)
            {
              dec_picture->ref_idx[LIST_0][img->block_x + i][img->block_y + j] = refframe;
            }
          
        }
      }
    }
  }
  else	//++ 计算非Direc类型块分割的前向参考帧（单帧参考，直接赋值）
  {
    for (j0=0; j0<4; j0+=step_v0)
    {
      for (i0=0; i0<4; i0+=step_h0)
      {
        k=2*(j0/2)+(i0/2);
        if ((currMB->b8pdir[k]==0 || currMB->b8pdir[k]==2) && currMB->b8mode[k]!=0)
        {
          for (j=j0; j<j0+step_v0;j++)
            for (i=i0; i<i0+step_h0;i++)
            {
              dec_picture->ref_idx[LIST_0][img->block_x + i][img->block_y + j] = 0;
            }
        }
      }
    }
  }
  
  //  If backward multiple ref. frames, read backward reference frame for the MB *********************************
  if(img->num_ref_idx_l1_active>1)	//++ 计算非Direc类型块分割的后向参考帧（多帧参考）
  {
    flag_mode = ( img->num_ref_idx_l1_active == 2 ? 1 : 0);	//++ 如果参考帧数只有两帧，则采用一比特定长编码。如果参考帧数大于2，则采用变长编码（哥伦布编码）
															//++ 这是因为如果参考帧数为2，采用FLC只需要一个比特就能表示参考帧号，而VLC可能需要三个比特（当参考帧号不为0时）
    currSE.type = SE_REFFRAME;
    dP = &(currSlice->partArr[partMap[SE_REFFRAME]]);
    if (active_pps->entropy_coding_mode_flag == UVLC || dP->bitstream->ei_flag)   currSE.mapping = linfo_ue;
    else                                                      currSE.reading = readRefFrame_CABAC;
    
    for (j0=0; j0<4; j0+=step_v0)
    {
      for (i0=0; i0<4; i0+=step_h0)
      {
        k=2*(j0/2)+(i0/2);
        if ((currMB->b8pdir[k]==1 || currMB->b8pdir[k]==2) && currMB->b8mode[k]!=0)
        {
          TRACE_STRING("ref_idx_l1");

          img->subblock_x = i0;
          img->subblock_y = j0;
          
          currSE.context = BType2CtxRef (currMB->b8mode[k]);
          if( (active_pps->entropy_coding_mode_flag == UVLC || dP->bitstream->ei_flag) && flag_mode )
          {
            currSE.len = 1;
            readSyntaxElement_FLC(&currSE, dP->bitstream);
            currSE.value1 = 1-currSE.value1;
          }
          else
          {
            currSE.value2 = LIST_1;
            dP->readSyntaxElement (&currSE,img,inp,dP);
          }
          refframe = currSE.value1;

          for (j=j0; j<j0+step_v0;j++)
          {
            for (i=i0; i<i0+step_h0;i++)
            {
              dec_picture->ref_idx[LIST_1][img->block_x + i][img->block_y + j] = refframe;
            }
          }
        }
      }
    }
  }
  else	//++ 计算非Direc类型块分割的后向参考帧（单帧参考，直接赋值）
  {
    for (j0=0; j0<4; j0+=step_v0)
    {
      for (i0=0; i0<4; i0+=step_h0)
      {
        k=2*(j0/2)+(i0/2);
        if ((currMB->b8pdir[k]==1 || currMB->b8pdir[k]==2) && currMB->b8mode[k]!=0)
        {
          for (j=j0; j<j0+step_v0;j++)
            for (i=i0; i<i0+step_h0;i++)
            {
              dec_picture->ref_idx[LIST_1][img->block_x + i][img->block_y + j] = 0;
            }
        }
      }
    }
  }

  //=====  READ FORWARD MOTION VECTORS =====
  currSE.type = SE_MVD;
  dP = &(currSlice->partArr[partMap[SE_MVD]]);

  if (active_pps->entropy_coding_mode_flag == UVLC || dP->bitstream->ei_flag) currSE.mapping = linfo_se;
  else                                                  currSE.reading = readMVD_CABAC;

  for (j0=0; j0<4; j0+=step_v0)
    for (i0=0; i0<4; i0+=step_h0)
    {
      k=2*(j0/2)+(i0/2);

      if ((currMB->b8pdir[k]==0 || currMB->b8pdir[k]==2) && (currMB->b8mode[k] !=0))//has forward vector
      {
        mv_mode  = currMB->b8mode[k];
        step_h   = BLOCK_STEP [mv_mode][0];
        step_v   = BLOCK_STEP [mv_mode][1];
        
        refframe = dec_picture->ref_idx[LIST_0][img->block_x+i0][img->block_y+j0];
        
        for (j=j0; j<j0+step_v0; j+=step_v)
        {
          for (i=i0; i<i0+step_h0; i+=step_h)
          {
            j4 = img->block_y+j;
            i4 = img->block_x+i;
            
            // first make mv-prediction
            SetMotionVectorPredictor (img, pmv, pmv+1, refframe, LIST_0, dec_picture->ref_idx, dec_picture->mv, i, j, 4*step_h, 4*step_v);

            for (k=0; k < 2; k++) 
            {
              TRACE_STRING("mvd_l0");

              img->subblock_x = i; // position used for context determination
              img->subblock_y = j; // position used for context determination
              currSE.value2 = k<<1; // identifies the component; only used for context determination
              dP->readSyntaxElement(&currSE,img,inp,dP);
              curr_mvd = currSE.value1; 
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
              if(en_ps.MVD)
              curr_mvd=MVD_decryption(curr_mvd,1);
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
              vec=curr_mvd+pmv[k];           /* find motion vector */
              
              for(ii=0;ii<step_h;ii++)
              {
                for(jj=0;jj<step_v;jj++)
                {
                  dec_picture->mv  [LIST_0][i4+ii][j4+jj][k] = vec;
                  currMB->mvd      [LIST_0][j+jj] [i+ii] [k] = curr_mvd;
                }
              }
            }
          }
        }
      }
      else if (currMB->b8mode[k=2*(j0/2)+(i0/2)]==0)      
      {  
        if (!img->direct_type)	//++ 时间直接模式
        {
          int list_offset = ((img->MbaffFrameFlag)&&(currMB->mb_field))? img->current_mb_nr%2 ? 4 : 2 : 0;
          int imgblock_y= ((img->MbaffFrameFlag)&&(currMB->mb_field))? (img->current_mb_nr%2) ? (img->block_y-4)/2:img->block_y/2 : img->block_y;

          int refList = (co_located_ref_idx[LIST_0][img->block_x+i0][imgblock_y+j0]== -1 ? LIST_1 : LIST_0);	//++ 如果“相对位置块”的前向参考帧不可用则采用其后向参考帧
          int ref_idx =  co_located_ref_idx[refList][img->block_x+i0][imgblock_y+j0];	//++ 此变量即为标准中的的refIdxCol
          
          if (ref_idx==-1)	//++ 该if语句对应标准中的公式8-188
          {
            for (j=j0; j<j0+step_v0; j++)
              for (i=i0; i<i0+step_h0; i++)
              {            
                dec_picture->ref_idx [LIST_1][img->block_x+i][img->block_y+j]=0;
                dec_picture->ref_idx [LIST_0][img->block_x+i][img->block_y+j]=0; 
                j4 = img->block_y+j;
                i4 = img->block_x+i;            
                for (ii=0; ii < 2; ii++) 
                {                                    
                  dec_picture->mv [LIST_0][i4][j4][ii]=0;
                  dec_picture->mv [LIST_1][i4][j4][ii]=0;                  
                }
              }
          }
          else 
          {        
            int mapped_idx=-1, iref;                             
            int j6;
                        
            for (iref=0;iref<min(img->num_ref_idx_l0_active,listXsize[LIST_0 + list_offset]);iref++)
            {
              
              if (dec_picture->ref_pic_num[img->current_slice_nr][LIST_0 + list_offset][iref]==co_located_ref_id[refList][img->block_x+i0][imgblock_y+j0])
              {
                mapped_idx=iref;
                break;
              }
              else //! invalid index. Default to zero even though this case should not happen
                mapped_idx=INVALIDINDEX;
            }
            
            if (INVALIDINDEX == mapped_idx)
            {
              error("temporal direct error\ncolocated block has ref that is unavailable",-1111);
            }
            
            
            for (j=j0; j<j0+step_v0; j++)
              for (i=i0; i<i0+step_h0; i++)
              {
                  mv_scale = img->mvscale[LIST_0 + list_offset][mapped_idx];

                  dec_picture->ref_idx [LIST_0][img->block_x+i][img->block_y+j] = mapped_idx;
                  dec_picture->ref_idx [LIST_1][img->block_x+i][img->block_y+j] = 0;
                  
                  j4 = img->block_y+j;
                  j6 = imgblock_y+j;
                  i4 = img->block_x+i;

                  for (ii=0; ii < 2; ii++) 
                  {              
                    //if (iTRp==0)
                    if (mv_scale == 9999 || listX[LIST_0+list_offset][mapped_idx]->is_long_term)
//                    if (mv_scale==9999 || Co_located->is_long_term)
                    {                      
                      dec_picture->mv  [LIST_0][i4][j4][ii]=co_located_mv[refList][i4][j6][ii];
                      dec_picture->mv  [LIST_1][i4][j4][ii]=0;
                    }
                    else
                    {
                      dec_picture->mv  [LIST_0][i4][j4][ii]=(mv_scale * co_located_mv[refList][i4][j6][ii] + 128 ) >> 8;
                      dec_picture->mv  [LIST_1][i4][j4][ii]=dec_picture->mv[LIST_0][i4][j4][ii] - co_located_mv[refList][i4][j6][ii];
                    }
                  }
              }
          }  
      } 
    }
  }
  
  //=====  READ BACKWARD MOTION VECTORS =====
  currSE.type = SE_MVD;
  dP          = &(currSlice->partArr[partMap[SE_MVD]]);

  if (active_pps->entropy_coding_mode_flag == UVLC || dP->bitstream->ei_flag) currSE.mapping = linfo_se;
  else                                                    currSE.reading = readMVD_CABAC;

  for (j0=0; j0<4; j0+=step_v0)
  {
    for (i0=0; i0<4; i0+=step_h0)
    {
      k=2*(j0/2)+(i0/2);
      if ((currMB->b8pdir[k]==1 || currMB->b8pdir[k]==2) && (currMB->b8mode[k]!=0))//has backward vector
      {
        mv_mode  = currMB->b8mode[k];
        step_h   = BLOCK_STEP [mv_mode][0];
        step_v   = BLOCK_STEP [mv_mode][1];
        
        refframe = dec_picture->ref_idx[LIST_1][img->block_x+i0][img->block_y+j0];

        for (j=j0; j<j0+step_v0; j+=step_v)
        {
          for (i=i0; i<i0+step_h0; i+=step_h)
          {
            j4 = img->block_y+j;
            i4 = img->block_x+i;
            
            // first make mv-prediction
            SetMotionVectorPredictor (img, pmv, pmv+1, refframe, LIST_1, dec_picture->ref_idx, dec_picture->mv, i, j, 4*step_h, 4*step_v);
            
            for (k=0; k < 2; k++) 
            {
              TRACE_STRING("mvd_l1");
              
              img->subblock_x = i; // position used for context determination
              img->subblock_y = j; // position used for context determination
              currSE.value2   = (k<<1) +1; // identifies the component; only used for context determination
              dP->readSyntaxElement(&currSE,img,inp,dP);
              curr_mvd = currSE.value1; 
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
              if(en_ps.MVD)
              curr_mvd=MVD_decryption(curr_mvd,1);
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
              vec=curr_mvd+pmv[k];           /* find motion vector */

              for(ii=0;ii<step_h;ii++)
              {
                for(jj=0;jj<step_v;jj++)
                {
                  dec_picture->mv  [LIST_1][i4+ii][j4+jj][k] = vec;
                  currMB->mvd      [LIST_1][j+jj] [i+ii] [k] = curr_mvd;
				}
              }
            }
          }
        }
      }
    }
  }
  // record reference picture Ids for deblocking decisions
 
  for(i4=img->block_x;i4<(img->block_x+4);i4++)
  for(j4=img->block_y;j4<(img->block_y+4);j4++)
  {
    if(dec_picture->ref_idx[LIST_0][i4][j4]>=0)
       dec_picture->ref_pic_id[LIST_0][i4][j4] = dec_picture->ref_pic_num[img->current_slice_nr][LIST_0 + list_offset][dec_picture->ref_idx[LIST_0][i4][j4]];
    else
       dec_picture->ref_pic_id[LIST_0][i4][j4] = INT64_MIN;
    if(dec_picture->ref_idx[LIST_1][i4][j4]>=0)
       dec_picture->ref_pic_id[LIST_1][i4][j4] = dec_picture->ref_pic_num[img->current_slice_nr][LIST_1 + list_offset][dec_picture->ref_idx[LIST_1][i4][j4]];  
    else
       dec_picture->ref_pic_id[LIST_1][i4][j4] = INT64_MIN;  
  }
}



/*!
 ************************************************************************
 * \brief
 *    Get the Prediction from the Neighboring Blocks for Number of Nonzero Coefficients 
 *    
 *    Luma Blocks
 ************************************************************************
 */
int predict_nnz(struct img_par *img, int i,int j)
{
  PixelPos pix;

  int pred_nnz = 0;
  int cnt      = 0;
  int mb_nr    = img->current_mb_nr;

  // left block
  getLuma4x4Neighbour(mb_nr, i, j, -1, 0, &pix);	//++ 当前块左边的块
/* to be inserted only for dp
  if (pix.available && active_pps->constrained_intra_pred_flag)
  {
    pix.available &= img->intra_block[pix.mb_addr];
  }
*/  
  if (pix.available)
  {
    pred_nnz = img->nz_coeff [pix.mb_addr ][pix.x][pix.y];
    cnt++;
  }

  // top block
  getLuma4x4Neighbour(mb_nr, i, j, 0, -1, &pix);	//++ 当前块上边的块
/* to be inserted only for dp
  if (pix.available && active_pps->constrained_intra_pred_flag)
  {
    pix.available &= img->intra_block[pix.mb_addr];
  }
*/  
  if (pix.available)
  {
    pred_nnz += img->nz_coeff [pix.mb_addr ][pix.x][pix.y];
    cnt++;
  }

  if (cnt==2)
  {
    pred_nnz++;
    pred_nnz/=cnt; 
  }

  return pred_nnz;
}


/*!
 ************************************************************************
 * \brief
 *    Get the Prediction from the Neighboring Blocks for Number of Nonzero Coefficients 
 *    
 *    Chroma Blocks   
 ************************************************************************
 */
int predict_nnz_chroma(struct img_par *img, int i,int j)
{
  PixelPos pix;

  int pred_nnz = 0;
  int cnt      =0;
  int mb_nr    = img->current_mb_nr;

  // left block
  getChroma4x4Neighbour(mb_nr, i%2, j-4, -1, 0, &pix);
/*  to be inserted only for dp
  if (pix.available && active_pps->constrained_intra_pred_flag)
  {
    pix.available &= img->intra_block[pix.mb_addr];
  }
*/  
  if (pix.available)
  {
    pred_nnz = img->nz_coeff [pix.mb_addr ][2 * (i/2) + pix.x][4 + pix.y];
    cnt++;
  }
  
  // top block
  getChroma4x4Neighbour(mb_nr, i%2, j-4, 0, -1, &pix);
/*  to be inserted only for dp
  if (pix.available && active_pps->constrained_intra_pred_flag)
  {
    pix.available &= img->intra_block[pix.mb_addr];
  }
*/  
  if (pix.available)
  {
    pred_nnz += img->nz_coeff [pix.mb_addr ][2 * (i/2) + pix.x][4 + pix.y];
    cnt++;
  }

  if (cnt==2)
  {
    pred_nnz++;
    pred_nnz/=cnt; 
  }

  return pred_nnz;
}


/*!
 ************************************************************************
 * \brief
 *    Reads coeff of an 4x4 block (CAVLC)
 *
 * \author
 *    Karl Lillevold <karll@real.com>
 *    contributions by James Au <james@ubvideo.com>
 ************************************************************************
 */

//++ 参见标准9.2.1
//++ levarr[k]为4*4块中按z-z顺序的第k个非零变换系数的值
//++ runarr[k]为4*4块中按z-z顺序的第k个变换系数前零的个数
void readCoeff4x4_CAVLC (struct img_par *img,struct inp_par *inp,
                        int block_type, 
                        int i, int j, int levarr[16], int runarr[16],
                        int *number_coefficients)
{
  int mb_nr = img->current_mb_nr;
  Macroblock *currMB = &img->mb_data[mb_nr];
  SyntaxElement currSE;
  Slice *currSlice = img->currentSlice;
  DataPartition *dP;
  int *partMap = assignSE2partition[currSlice->dp_mode];


  int k, code, vlcnum;
  int numcoeff, numtrailingones, numcoeff_vlc;
  int level_two_or_higher;
  int numones, totzeros, level, cdc=0, cac=0;
  int zerosleft, ntr, dptype = 0;
  int max_coeff_num = 0, nnz;	//++ 此处的nnz就相当于9.2.1节中的nC
  char type[15];
  int incVlc[] = {0,3,6,12,24,48,32768};    // maximum vlc = 6
											//++ 对应的部分位于标准的"7.3.5.3.1 Residual block CAVLC syntax"，
											//++ 它表明当前编解码的量化级(Level)的大小超过incVlc[suffixLength]时，
  											//++ suffixLength应该加1。incVlc[0] = 0，表示编解码一个量化级之后，
  											//++ suffixLength至少会增加一，因为量化级的大小肯定大于0；
  											//++ incVlc[6] = 32768表示suffixLength最大只能为6，因为量化级的大
  											//++ 小不会超过2^15，因而suffixLength只能维持为6；
  											//++ incVlc[i] = 3 << (i - 1)，i = 1…5。所以它的作用就是控制
  											//++ suffixLength的增长速度。
  numcoeff = 0;

  switch (block_type)
  {
  case LUMA:
    max_coeff_num = 16;
    sprintf(type, "%s", "Luma");
    if (IS_INTRA (currMB))
    {
      dptype = SE_LUM_AC_INTRA;
    }
    else
    {
      dptype = SE_LUM_AC_INTER;
    }
    break;
  case LUMA_INTRA16x16DC:
    max_coeff_num = 16;
    sprintf(type, "%s", "Lum16DC");
    dptype = SE_LUM_DC_INTRA;
    break;
  case LUMA_INTRA16x16AC:
    max_coeff_num = 15;
    sprintf(type, "%s", "Lum16AC");
    dptype = SE_LUM_AC_INTRA;
    break;

  case CHROMA_DC:
    max_coeff_num = 4;
    cdc = 1;

    sprintf(type, "%s", "ChrDC");
    if (IS_INTRA (currMB))
    {
      dptype = SE_CHR_DC_INTRA;
    }
    else
    {
      dptype = SE_CHR_DC_INTER;
    }
    break;
  case CHROMA_AC:
    max_coeff_num = 15;
    cac = 1;
    sprintf(type, "%s", "ChrAC");
    if (IS_INTRA (currMB))
    {
      dptype = SE_CHR_AC_INTRA;
    }
    else
    {
      dptype = SE_CHR_AC_INTER;
    }
    break;
  default:
    error ("readCoeff4x4_CAVLC: invalid block type", 600);
    break;
  }

  currSE.type = dptype;
  dP = &(currSlice->partArr[partMap[dptype]]);

  img->nz_coeff[img->current_mb_nr][i][j] = 0;


  if (!cdc)	//++ 判断色度直流分量
  {
    // luma or chroma AC
    if (!cac)	//++ 判断色度交流分量
    {
      nnz = predict_nnz(img, i, j);	//++ 得到非零系数个数的预测值
    }
    else
    {
      nnz = predict_nnz_chroma(img, i, j);
    }

    if (nnz < 2)
    {
      numcoeff_vlc = 0;	//++ 即使用表9-5第三列（变长码表）
    }
    else if (nnz < 4)
    {
      numcoeff_vlc = 1;	//++ 即使用表9-5第四列（变长码表）
    }
    else if (nnz < 8)
    {
      numcoeff_vlc = 2;	//++ 即使用表9-5第五列（变长码表）
    }
    else //
    {
      numcoeff_vlc = 3;	//++ 即使用表9-5第六列（6比特定长码表）
    }

    currSE.value1 = numcoeff_vlc;

    readSyntaxElement_NumCoeffTrailingOnes(&currSE, dP, type);

    numcoeff =  currSE.value1;
    numtrailingones =  currSE.value2;

    img->nz_coeff[img->current_mb_nr][i][j] = numcoeff;
  }
  else
  {
    // chroma DC
    readSyntaxElement_NumCoeffTrailingOnesChromaDC(&currSE, dP);

    numcoeff =  currSE.value1;
    numtrailingones =  currSE.value2;
  }


  for (k = 0; k < max_coeff_num; k++)
  {
    levarr[k] = 0;
    runarr[k] = 0;
  }

  numones = numtrailingones;
  *number_coefficients = numcoeff;

  if (numcoeff)
  {
    if (numtrailingones)
    {

      currSE.len = numtrailingones;

#if TRACE
      snprintf(currSE.tracestring, 
        TRACESTRING_SIZE, "%s trailing ones sign (%d,%d)", type, i, j);
#endif

      readSyntaxElement_FLC (&currSE, dP->bitstream);

      code = currSE.inf;
      ntr = numtrailingones;
      for (k = numcoeff-1; k > numcoeff-1-numtrailingones; k--)
      {
        ntr --;
        if ((code>>ntr)&1)
          levarr[k] = -1;
        else
          levarr[k] = 1;
      }
    }

    // decode levels
    level_two_or_higher = 1;
    if (numcoeff > 3 && numtrailingones == 3)
      level_two_or_higher = 0;

      if (numcoeff > 10 && numtrailingones < 3)
          vlcnum = 1;
      else
          vlcnum = 0;

    for (k = numcoeff - 1 - numtrailingones; k >= 0; k--)
    {

#if TRACE
      snprintf(currSE.tracestring, 
        TRACESTRING_SIZE, "%s lev (%d,%d) k=%d vlc=%d ", type,
          i, j, k, vlcnum);
#endif

      if (vlcnum == 0)
          readSyntaxElement_Level_VLC0(&currSE, dP);
      else
          readSyntaxElement_Level_VLCN(&currSE, vlcnum, dP);

      if (level_two_or_higher)
      {
          if (currSE.inf > 0)
          currSE.inf ++;
          else
          currSE.inf --;
          level_two_or_higher = 0;
      }

      level = levarr[k] = currSE.inf;
      if (abs(level) == 1)
        numones ++;

      // update VLC table
      if (abs(level)>incVlc[vlcnum])
        vlcnum++;

      if (k == numcoeff - 1 - numtrailingones && abs(level)>3)
        vlcnum = 2;

    }
    
    if (numcoeff < max_coeff_num)
    {
      // decode total run
      vlcnum = numcoeff-1;
      currSE.value1 = vlcnum;

#if TRACE
      snprintf(currSE.tracestring, 
        TRACESTRING_SIZE, "%s totalrun (%d,%d) vlc=%d ", type, i,j, vlcnum);
#endif
      if (cdc)
        readSyntaxElement_TotalZerosChromaDC(&currSE, dP);
      else
        readSyntaxElement_TotalZeros(&currSE, dP);

      totzeros = currSE.value1;
    }
    else
    {
      totzeros = 0;
    }

    // decode run before each coefficient
    zerosleft = totzeros;
    i = numcoeff-1;
    if (zerosleft > 0 && i > 0)
    {
      do 
      {
        // select VLC for runbefore
        vlcnum = zerosleft - 1;
        if (vlcnum > RUNBEFORE_NUM-1)
          vlcnum = RUNBEFORE_NUM-1;

        currSE.value1 = vlcnum;
#if TRACE
        snprintf(currSE.tracestring, 
          TRACESTRING_SIZE, "%s run (%d,%d) k=%d vlc=%d ",
            type, i, j, i, vlcnum);
#endif

        readSyntaxElement_Run(&currSE, dP);
        runarr[i] = currSE.value1;

        zerosleft -= runarr[i];
        i --;
      } while (zerosleft != 0 && i != 0);
    }
    runarr[i] = zerosleft;

  } // if numcoeff
}



/*!
 ************************************************************************
 * \brief
 *    Get coded block pattern and coefficients (run/level)
 *    from the NAL
 ************************************************************************
 */
void readCBPandCoeffsFromNAL(struct img_par *img,struct inp_par *inp)
{
  static int count=0;///////////////////////////////////////////
  int i,j,k;
  int level;
  int mb_nr = img->current_mb_nr;
  int ii,jj;
  int i1,j1, m2,jg2;
  Macroblock *currMB = &img->mb_data[mb_nr];
  int cbp;
  SyntaxElement currSE;
  Slice *currSlice = img->currentSlice;
  DataPartition *dP;
  int *partMap = assignSE2partition[currSlice->dp_mode];
  int iii,jjj;
  int coef_ctr, i0, j0, b8;
  int ll;
  int block_x,block_y;
  int start_scan;
  int uv;
  int qp_uv;
  int run, len;
  int levarr[16], runarr[16], numcoeff;

  int qp_per    = (img->qp-MIN_QP)/6;
  int qp_rem    = (img->qp-MIN_QP)%6;
  int qp_per_uv = QP_SCALE_CR[img->qp-MIN_QP]/6;
  int qp_rem_uv = QP_SCALE_CR[img->qp-MIN_QP]%6;
  int smb       = ((img->type==SP_SLICE) && IS_INTER (currMB)) || (img->type == SI_SLICE && currMB->mb_type == SI4MB);

  // QPI
  qp_uv = img->qp + active_pps->chroma_qp_index_offset;
  qp_uv = Clip3(0, 51, qp_uv);
  qp_per_uv = QP_SCALE_CR[qp_uv-MIN_QP]/6;
  qp_rem_uv = QP_SCALE_CR[qp_uv-MIN_QP]%6;

  // read CBP if not new intra mode
  if (!IS_NEWINTRA (currMB))	//++ 如果当前宏块不是帧内16*16或者I_PCM类型宏块则从码流中读取CBP
  {
    if (IS_OLDINTRA (currMB) || currMB->mb_type == SI4MB )   currSE.type = SE_CBP_INTRA;
    else                        currSE.type = SE_CBP_INTER;

    dP = &(currSlice->partArr[partMap[currSE.type]]);
    
    if (active_pps->entropy_coding_mode_flag == UVLC || dP->bitstream->ei_flag)
    {
      if (IS_OLDINTRA (currMB) || currMB->mb_type == SI4MB)  currSE.mapping = linfo_cbp_intra;
      else                       currSE.mapping = linfo_cbp_inter;
    }
    else
    {
      currSE.reading = readCBP_CABAC;
    }

    TRACE_STRING("coded_block_pattern");

    dP->readSyntaxElement(&currSE,img,inp,dP);
    currMB->cbp = cbp = currSE.value1;
    // Delta quant only if nonzero coeffs
    if (cbp !=0)
    {
      if (IS_INTER (currMB))  currSE.type = SE_DELTA_QUANT_INTER;
      else                    currSE.type = SE_DELTA_QUANT_INTRA;

      dP = &(currSlice->partArr[partMap[currSE.type]]);
      
      if (active_pps->entropy_coding_mode_flag == UVLC || dP->bitstream->ei_flag)
      {
        currSE.mapping = linfo_se;
      }
      else
                currSE.reading= readDquant_CABAC; //gabi

      TRACE_STRING("mb_qp_delta");

      dP->readSyntaxElement(&currSE,img,inp,dP);
      currMB->delta_quant = currSE.value1;
      img->qp= (img->qp-MIN_QP+currMB->delta_quant+(MAX_QP-MIN_QP+1))%(MAX_QP-MIN_QP+1)+MIN_QP;
    }
  }
  else
  {
    cbp = currMB->cbp;
  }

  for (i=0;i<BLOCK_SIZE;i++)
    for (j=0;j<BLOCK_SIZE;j++)
      for(iii=0;iii<BLOCK_SIZE;iii++)
        for(jjj=0;jjj<BLOCK_SIZE;jjj++)
          img->cof[i][j][iii][jjj]=0;// reset luma coeffs


  if (IS_NEWINTRA (currMB)) // read DC coeffs for new intra modes	//++ 当前宏块是帧内16*16或者I_PCM类型宏块
  {
    currSE.type = SE_DELTA_QUANT_INTRA;

    dP = &(currSlice->partArr[partMap[currSE.type]]);
    
    if (active_pps->entropy_coding_mode_flag == UVLC || dP->bitstream->ei_flag)
    {
      currSE.mapping = linfo_se;
    }
    else
    {
      currSE.reading= readDquant_CABAC;
    }
#if TRACE
    snprintf(currSE.tracestring, TRACESTRING_SIZE, "Delta quant ");
#endif
    dP->readSyntaxElement(&currSE,img,inp,dP);
    currMB->delta_quant = currSE.value1;
    img->qp= (img->qp-MIN_QP+currMB->delta_quant+(MAX_QP-MIN_QP+1))%(MAX_QP-MIN_QP+1)+MIN_QP;

    for (i=0;i<BLOCK_SIZE;i++)
      for (j=0;j<BLOCK_SIZE;j++)
        img->ipredmode[img->block_x+i][img->block_y+j]=DC_PRED;

    if (active_pps->entropy_coding_mode_flag == UVLC)
    {
      readCoeff4x4_CAVLC(img, inp, LUMA_INTRA16x16DC, 0, 0,
                          levarr, runarr, &numcoeff);
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
     if(en_ps.DCT)
      DCT_no_coeff_decryption(levarr,runarr,0,en_ps.DCT);
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      coef_ctr=-1;
      level = 1;                            // just to get inside the loop
      for(k = 0; k < numcoeff; k++)
      {
        if (levarr[k] != 0)                     // leave if len=1
        {
          coef_ctr=coef_ctr+runarr[k]+1;

          if ((img->structure == FRAME) && (!currMB->mb_field)) 
          {
            i0=SNGL_SCAN[coef_ctr][0];
            j0=SNGL_SCAN[coef_ctr][1];
          }
          else 
          { // Alternate scan for field coding
            i0=FIELD_SCAN[coef_ctr][0];
            j0=FIELD_SCAN[coef_ctr][1];
          }

          img->cof[i0][j0][0][0]=levarr[k];// add new intra DC coeff
        }
      }
    }
    else
    {

      currSE.type = SE_LUM_DC_INTRA;
      dP = &(currSlice->partArr[partMap[currSE.type]]);

      currSE.context      = LUMA_16DC;
      currSE.type         = SE_LUM_DC_INTRA;
      img->is_intra_block = 1;

      if (active_pps->entropy_coding_mode_flag == UVLC || dP->bitstream->ei_flag)
      {
        currSE.mapping = linfo_levrun_inter;
      }
      else
      {
        currSE.reading = readRunLevel_CABAC;
      }



      coef_ctr=-1;
      level = 1;                            // just to get inside the loop
      for(k=0;(k<17) && (level!=0);k++)
      {
#if TRACE
        snprintf(currSE.tracestring, TRACESTRING_SIZE, "DC luma 16x16 ");
#endif
        dP->readSyntaxElement(&currSE,img,inp,dP);
        level = currSE.value1;
        run   = currSE.value2;
        len   = currSE.len;

        if (level != 0)                     // leave if len=1
        {
          coef_ctr=coef_ctr+run+1;

          if ((img->structure == FRAME) && (!currMB->mb_field)) 
          {
            i0=SNGL_SCAN[coef_ctr][0];
            j0=SNGL_SCAN[coef_ctr][1];
          }
          else 
          { // Alternate scan for field coding
            i0=FIELD_SCAN[coef_ctr][0];
            j0=FIELD_SCAN[coef_ctr][1];
          }

          img->cof[i0][j0][0][0]=level;// add new intra DC coeff
        }
      }
    }

    itrans_2(img);// transform new intra DC	//++ 亮度DC系数反Hadamard变换
  }

  qp_per    = (img->qp-MIN_QP)/6;
  qp_rem    = (img->qp-MIN_QP)%6;
  qp_uv = img->qp + active_pps->chroma_qp_index_offset;
  qp_uv = Clip3(0, 51, qp_uv);
  qp_per_uv = QP_SCALE_CR[qp_uv-MIN_QP]/6;
  qp_rem_uv = QP_SCALE_CR[qp_uv-MIN_QP]%6;
  currMB->qp = img->qp;

  // luma coefficients
  //++ ************* 读取一个宏块中每个4*4块亮度系数 *************
  //++ 4*4块的处理顺序为：
  //++		***************
  //++ 		* 1  2  5  6  *
  //++ 		* 3  4  7  8  *
  //++ 		* 9  10 13 14 *
  //++ 		* 11 12 15 16 *
  //++		***************
  //++ 4*4块中变换系数的读取顺序为：
  //++		***************
  //++ 		* 1  2  6  7  *
  //++ 		* 3  5  8  13 *
  //++ 		* 4  9  12 14 *
  //++ 		* 10 11 15 16 *
  //++		***************
  for (block_y=0; block_y < 4; block_y += 2) /* all modes */	//++++++++++++++++++++++++++++++++++++++++
  {																//++   该双重循环为16*16块中的4个8*8块
    for (block_x=0; block_x < 4; block_x += 2)					//++++++++++++++++++++++++++++++++++++++++
    {

      b8 = 2*(block_y/2) + block_x/2;
      if (active_pps->entropy_coding_mode_flag == UVLC)
      {
        for (j=block_y; j < block_y+2; j++)		//++++++++++++++++++++++++++++++++++++++
        {										//++   该双重循环为8*8块中的4个4*4块
          for (i=block_x; i < block_x+2; i++)	//++++++++++++++++++++++++++++++++++++++
          {
            ii = block_x/2; jj = block_y/2;
            b8 = 2*jj+ii;

            if (cbp & (1<<b8))  /* are there any coeff in current block at all */
            {
              if (IS_NEWINTRA(currMB))
              {
				//++ 读取4*4亮度块的变换系数，按z-z顺序将第k个非零变换系数的值放入levarr[k]中
				//++ 将第k个变换系数前零的个数放入runarr[k]中
                readCoeff4x4_CAVLC(img, inp, LUMA_INTRA16x16AC, i, j,
                                    levarr, runarr, &numcoeff);
                start_scan = 1;
              }
              else
              {
				//++ 读取4*4亮度块的变换系数，按z-z顺序将第k个非零变换系数的值放入levarr[k]中
				//++ 将第k个变换系数前零的个数放入runarr[k]中
                readCoeff4x4_CAVLC(img, inp, LUMA, i, j,
                                    levarr, runarr, &numcoeff);
                start_scan = 0;
              }
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				 if(en_ps.DCT)
				  DCT_no_coeff_decryption(levarr,runarr,0,en_ps.DCT);
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
              coef_ctr = start_scan-1;	//++ coef_ctr为系数计数，即指明当前处理的是Z-Z顺序中的第几个变换系数
              for (k = 0; k < numcoeff; k++)
              {
                if (levarr[k] != 0)
                {
                  coef_ctr += runarr[k]+1;

                  if ((img->structure == FRAME) && (!currMB->mb_field)) 
                  {
                    i0=SNGL_SCAN[coef_ctr][0];
                    j0=SNGL_SCAN[coef_ctr][1];
                  }
                  else { // Alternate scan for field coding
                    i0=FIELD_SCAN[coef_ctr][0];
                    j0=FIELD_SCAN[coef_ctr][1];
                  }
                  currMB->cbp_blk      |= 1 << ((j<<2) + i) ;	//++ 当前宏块中当前4*4块的cbp标志位的值

                  img->cof[i][j][i0][j0]= levarr[k]*dequant_coef[qp_rem][i0][j0]<<qp_per;	//++ i,j标识16*16块中的哪一个4*4块;i0,j0标识4*4块内的哪一个位置的变换系数
                }
              }
            }
            else
            {
              img->nz_coeff[img->current_mb_nr][i][j] = 0;
            }
          }
        }
      } // VLC
      else
      { 
          b8 = 2*(block_y/2) + block_x/2;
          // CABAC 
          for (j=block_y; j < block_y+2; j++)
          {
            //jj=j/2;
            for (i=block_x; i < block_x+2; i++)
            {
              //ii = i/2;
              //b8 = 2*jj+ii;

              if (IS_NEWINTRA (currMB))   start_scan = 1; // skip DC coeff
              else                        start_scan = 0; // take all coeffs

              img->subblock_x = i; // position for coeff_count ctx
              img->subblock_y = j; // position for coeff_count ctx
              if (cbp & (1<<b8))  // are there any coeff in current block at all
              {
                coef_ctr = start_scan-1;
                level    = 1;
                for(k=start_scan;(k<17) && (level!=0);k++)
                {
                 /*
                  * make distinction between INTRA and INTER coded
                  * luminance coefficients
                  */
                  currSE.context      = (IS_NEWINTRA(currMB) ? LUMA_16AC : LUMA_4x4);
                  currSE.type         = (IS_INTRA(currMB) ?
                                        (k==0 ? SE_LUM_DC_INTRA : SE_LUM_AC_INTRA) :
                                        (k==0 ? SE_LUM_DC_INTER : SE_LUM_AC_INTER));
                  img->is_intra_block = IS_INTRA(currMB);
                  
#if TRACE
                  sprintf(currSE.tracestring, "Luma sng ");
#endif
                  dP = &(currSlice->partArr[partMap[currSE.type]]);
                  
                  if (active_pps->entropy_coding_mode_flag == UVLC || dP->bitstream->ei_flag)  currSE.mapping = linfo_levrun_inter;
                  else                                                     currSE.reading = readRunLevel_CABAC;
                  
                  dP->readSyntaxElement(&currSE,img,inp,dP);
                  level = currSE.value1;
                  run   = currSE.value2;
                  len   = currSE.len;
                  
                  if (level != 0)    /* leave if len=1 */
                  {
                    coef_ctr             += run+1;

                    if ((img->structure == FRAME) && (!currMB->mb_field)) 
                    {
                      i0=SNGL_SCAN[coef_ctr][0];
                      j0=SNGL_SCAN[coef_ctr][1];
                    }
                    else { // Alternate scan for field coding
                      i0=FIELD_SCAN[coef_ctr][0];
                      j0=FIELD_SCAN[coef_ctr][1];
                    }
                    currMB->cbp_blk      |= 1 << ((j<<2) + i) ;
                    img->cof[i][j][i0][j0]= level*dequant_coef[qp_rem][i0][j0]<<qp_per;
                  }
                }
              }
            }
          }
      } 
    }
  }

  for (j=4;j<6;j++) // reset all chroma coeffs before read
    for (i=0;i<4;i++)
      for (iii=0;iii<4;iii++)
        for (jjj=0;jjj<4;jjj++)
          img->cof[i][j][iii][jjj]=0;

  m2 =img->mb_x*2;
  jg2=img->mb_y*2;

  // chroma 2x2 DC coeff
  if(cbp>15)
  {
    for (ll=0;ll<3;ll+=2)
    {
      for (i=0;i<4;i++)
        img->cofu[i]=0;


      if (active_pps->entropy_coding_mode_flag == UVLC)
      {

        readCoeff4x4_CAVLC(img, inp, CHROMA_DC, 0, 0,
                            levarr, runarr, &numcoeff);
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	   if(en_ps.DCT)
		DCT_no_coeff_decryption(levarr,runarr,1,en_ps.DCT);
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        coef_ctr=-1;
        level=1;
        for(k = 0; k < numcoeff; k++)
        {
          if (levarr[k] != 0)
          {
            currMB->cbp_blk |= 0xf0000 << (ll<<1) ;
            coef_ctr=coef_ctr+runarr[k]+1;
            img->cofu[coef_ctr]=levarr[k];
          }
        }
      }
      else
      {
        coef_ctr=-1;
        level=1;
        for(k=0;(k<5)&&(level!=0);k++)
        {
          currSE.context      = CHROMA_DC;
          currSE.type         = (IS_INTRA(currMB) ? SE_CHR_DC_INTRA : SE_CHR_DC_INTER);
          img->is_intra_block =  IS_INTRA(currMB);
          img->is_v_block     = ll;

#if TRACE
          snprintf(currSE.tracestring, TRACESTRING_SIZE, " 2x2 DC Chroma ");
#endif
          dP = &(currSlice->partArr[partMap[currSE.type]]);
        
          if (active_pps->entropy_coding_mode_flag == UVLC || dP->bitstream->ei_flag)
            currSE.mapping = linfo_levrun_c2x2;
          else
            currSE.reading = readRunLevel_CABAC;

          dP->readSyntaxElement(&currSE,img,inp,dP);
          level = currSE.value1;
          run = currSE.value2;
          len = currSE.len;
          if (level != 0)
          {
            currMB->cbp_blk |= 0xf0000 << (ll<<1) ;
            coef_ctr=coef_ctr+run+1;
            // Bug: img->cofu has only 4 entries, hence coef_ctr MUST be <4 (which is
            // caught by the assert().  If it is bigger than 4, it starts patching the
            // img->predmode pointer, which leads to bugs later on.
            //
            // This assert() should be left in the code, because it captures a very likely
            // bug early when testing in error prone environments (or when testing NAL
            // functionality).
            assert (coef_ctr < 4);
            img->cofu[coef_ctr]=level;
          }
        }
      }

      if (smb) // check to see if MB type is SPred or SIntra4x4 
      {
         img->cof[0+ll][4][0][0]=img->cofu[0];   img->cof[1+ll][4][0][0]=img->cofu[1];
         img->cof[0+ll][5][0][0]=img->cofu[2];   img->cof[1+ll][5][0][0]=img->cofu[3];
      }
      else
      {
        for (i=0;i<4;i++)
          img->cofu[i]*=dequant_coef[qp_rem_uv][0][0]<<qp_per_uv;							//++ 色度DC系数反量化
        img->cof[0+ll][4][0][0]=(img->cofu[0]+img->cofu[1]+img->cofu[2]+img->cofu[3])>>1;	//++++++++++++++++++
        img->cof[1+ll][4][0][0]=(img->cofu[0]-img->cofu[1]+img->cofu[2]-img->cofu[3])>>1;	//++ 色度DC系数反
        img->cof[0+ll][5][0][0]=(img->cofu[0]+img->cofu[1]-img->cofu[2]-img->cofu[3])>>1;	//++ Hadamard变换
        img->cof[1+ll][5][0][0]=(img->cofu[0]-img->cofu[1]-img->cofu[2]+img->cofu[3])>>1;	//+++++++++++++++++
      }
    }
  }

  // chroma AC coeff, all zero from start_scan
  if (cbp<=31)
    for (j=4; j < 6; j++)
      for (i=0; i < 4; i++)
        img->nz_coeff [img->current_mb_nr ][i][j]=0;


  // chroma AC coeff, all zero from start_scan
  uv=-1;
  if (cbp>31)
  {
    block_y=4;
    for (block_x=0; block_x < 4; block_x += 2)
    {
      for (j=block_y; j < block_y+2; j++)
      {
        jj=j/2;
        j1=j-4;
        for (i=block_x; i < block_x+2; i++)
        {

          ii=i/2;
          i1=i%2;

          if (active_pps->entropy_coding_mode_flag == UVLC)
          {
            readCoeff4x4_CAVLC(img, inp, CHROMA_AC, i, j,
                                levarr, runarr, &numcoeff);
     /////////////////////////////////////////////////////////////////////
	      if(en_ps.DCT)
		  DCT_no_coeff_decryption(levarr,runarr,0,en_ps.DCT);
     /////////////////////////////////////////////////////////////////////
            coef_ctr=0;
            level=1;
            uv++;
            for(k = 0; k < numcoeff;k++)
            {
              if (levarr[k] != 0)
              {
                currMB->cbp_blk |= 1 << (16 + (j1<<1) + i1 + (block_x<<1) ) ;
                coef_ctr=coef_ctr+runarr[k]+1;

                if ((img->structure == FRAME) && (!currMB->mb_field)) 
                {
                  i0=SNGL_SCAN[coef_ctr][0];
                  j0=SNGL_SCAN[coef_ctr][1];
                }
                else { // Alternate scan for field coding
                  i0=FIELD_SCAN[coef_ctr][0];
                  j0=FIELD_SCAN[coef_ctr][1];
                }

                img->cof[i][j][i0][j0]=levarr[k]*dequant_coef[qp_rem_uv][i0][j0]<<qp_per_uv;
              }
            }
          }

          else
          {
            coef_ctr=0;
            level=1;
            uv++;

            img->subblock_y = j/5;
            img->subblock_x = i%2;

            for(k=0;(k<16)&&(level!=0);k++)
            {
              currSE.context      = CHROMA_AC;
              currSE.type         = (IS_INTRA(currMB) ? SE_CHR_AC_INTRA : SE_CHR_AC_INTER);
              img->is_intra_block =  IS_INTRA(currMB);
              img->is_v_block     = (uv>=4);

#if TRACE
              snprintf(currSE.tracestring, TRACESTRING_SIZE, " AC Chroma ");
#endif
              dP = &(currSlice->partArr[partMap[currSE.type]]);
            
              if (active_pps->entropy_coding_mode_flag == UVLC || dP->bitstream->ei_flag)
                currSE.mapping = linfo_levrun_inter;
              else
                currSE.reading = readRunLevel_CABAC;
                       
              dP->readSyntaxElement(&currSE,img,inp,dP);
              level = currSE.value1;
              run = currSE.value2;
              len = currSE.len;

              if (level != 0)
              {
                currMB->cbp_blk |= 1 << (16 + (j1<<1) + i1 + (block_x<<1) ) ;
                coef_ctr=coef_ctr+run+1;

                if ((img->structure == FRAME) && (!currMB->mb_field)) 
                {
                  i0=SNGL_SCAN[coef_ctr][0];
                  j0=SNGL_SCAN[coef_ctr][1];
                }
                else { // Alternate scan for field coding
                  i0=FIELD_SCAN[coef_ctr][0];
                  j0=FIELD_SCAN[coef_ctr][1];
                }
                img->cof[i][j][i0][j0]=level*dequant_coef[qp_rem_uv][i0][j0]<<qp_per_uv;
              }
            }
          }
        }
      }
    }
  }
}


/*!
 ************************************************************************
 * \brief
 *    Copy IPCM coefficients to decoded picture buffer and set parameters for this MB
 *    (for IPCM CABAC and IPCM CAVLC  28/11/2003)
 *
 * \author
 *    Dong Wang <Dong.Wang@bristol.ac.uk>
 ************************************************************************
 */

void decode_ipcm_mb(struct img_par *img)
{
  int i,j;

  Macroblock *currMb = &img->mb_data[img->current_mb_nr];    

  //Copy coefficents to decoded picture buffer
  //IPCM coefficents are stored in img->cof which is set in function readIPCMcoeffsFromNAL()

  for(i=0;i<16;i++)
    for(j=0;j<16;j++)
      dec_picture->imgY[img->pix_y+i][img->pix_x+j]=img->cof[i/4][j/4][i%4][j%4];

  for(i=0;i<8;i++)
    for(j=0;j<8;j++)
      dec_picture->imgUV[0][img->pix_c_y+i][img->pix_c_x+j]=img->cof[i/4][j/4+4][i%4][j%4];

  for(i=0;i<8;i++)
    for(j=0;j<8;j++)
      dec_picture->imgUV[1][img->pix_c_y+i][img->pix_c_x+j]=img->cof[i/4+2][j/4+4][i%4][j%4];



  //For Deblocking Filter  16/08/2003
  if (currMb->mb_type==IPCM)
    currMb->qp=0;

  //For CAVLC
  //Set the nz_coeff to 16. 
  //These parameters are to be used in CAVLC decoding of neighbour blocks
  for(i=0;i<4;i++)
    for (j=0;j<6;j++)
      img->nz_coeff[img->current_mb_nr][i][j]=16;


  //For CABAC decoding of MB skip flag 
  if (currMb->mb_type==IPCM)
    currMb->skip_flag=1;

  //for Loop filter CABAC
  if (currMb->mb_type==IPCM)
    currMb->cbp_blk=0xFFFF;

  //For CABAC decoding of Dquant
  last_dquant=0;
}

/*!
 ************************************************************************
 * \brief
 *    decode one macroblock
 ************************************************************************
 */

/***************************************************************************************************************
decode_one_slice()函数的流程（decode_one_macroblock()是decode_one_slice()一个比较重要的函数）：
decode_one_macroblock()首先针对亮度块解码，之后对色度块解码。首先判断宏块是否采用16x16帧内预测，如果是，
则通过函数intrapred_luma_16x16()解码。然后对一个16x16宏块中的每个4x4块分别解码（每个块的解码顺序见毕厚杰书
的图6.49（P118）），直到所有块解码完成为止。解码4x4块时，首先判断是否为帧内编码，如果是，则通过函数intrapred
对4x4块帧内解码。如果是帧间编码，则还要判断是P或是B帧。找到参考帧的匹配块后，用get_block()像素内插恢复像素值。
itrans()对当前块的残差进行反变换，然后再将反变换值和预测值相加得到重构值。接着，将该重构值赋给dec_picture->imagY
保存。其中de_picture是一个用来存放解码值的解构体，其成员imgy存放亮度解码值。最后是色度块解码。色度块解码和亮度块
解码流程基本相同。
******************************************************************************************************************/
int decode_one_macroblock(struct img_par *img,struct inp_par *inp)
{
  int tmp_block[BLOCK_SIZE][BLOCK_SIZE];
  int tmp_blockbw[BLOCK_SIZE][BLOCK_SIZE];
  int i=0,j=0,k,l,ii=0,jj=0,i1=0,j1=0,j4=0,i4=0;
  int uv, hv;
  int vec1_x=0,vec1_y=0,vec2_x=0,vec2_y=0;	//++ 当前4*4块的参考块左上角像素在本图像帧中的坐标（以1/4像素距离为单位）
  int ioff,joff;
  int block8x8;   // needed for ABT

  int bw_pred=0, fw_pred=0, pred, ifx;
  int ii0,jj0,ii1,jj1,if1,jf1,if0,jf0;
  int mv_mul,f1,f2,f3,f4;

  const byte decode_block_scan[16] = {0,1,4,5,2,3,6,7,8,9,12,13,10,11,14,15};	//++ 参见6.4.3	Inverse 4x4 luma block scanning process

  Macroblock *currMB   = &img->mb_data[img->current_mb_nr];
  int ref_idx, fw_refframe=-1, bw_refframe=-1, mv_mode, pred_dir, intra_prediction; // = currMB->ref_frame;
  int fw_ref_idx=-1, bw_ref_idx=-1;

  int  *** mv_array, ***fw_mv_array, ***bw_mv_array;

  int mv_scale;

  int mb_nr             = img->current_mb_nr;
  int smb       = ((img->type==SP_SLICE) && IS_INTER (currMB)) || (img->type == SI_SLICE && currMB->mb_type == SI4MB);
  int list_offset;
  int max_y_cr;

  StorablePicture **list;

  int jf;

  int fw_rFrame=-1,bw_rFrame=-1;
  int pmvfw[2]={0,0},pmvbw[2]={0,0};              

  int direct_pdir=-1;

  int curr_mb_field = ((img->MbaffFrameFlag)&&(currMB->mb_field));
  
  byte **     moving_block;
  int ****     co_located_mv;
  int ***     co_located_ref_idx;
  int64 ***    co_located_ref_id;

  if(currMB->mb_type==IPCM)
  {
    //copy readed data into imgY and set parameters
    decode_ipcm_mb(img);
    return 0;
  }

//////////////////////////

  // find out the correct list offsets
  if (curr_mb_field)
  {
    if(mb_nr%2)
    {
      list_offset = 4; // top field mb
      moving_block = Co_located->bottom_moving_block;
      co_located_mv = Co_located->bottom_mv;
      co_located_ref_idx = Co_located->bottom_ref_idx;
      co_located_ref_id = Co_located->bottom_ref_pic_id;
    }
    else
    {
      list_offset = 2; // bottom field mb
      moving_block = Co_located->top_moving_block;
      co_located_mv = Co_located->top_mv;
      co_located_ref_idx = Co_located->top_ref_idx;
      co_located_ref_id = Co_located->top_ref_pic_id;
    }
    max_y_cr = dec_picture->size_y_cr/2-1;
  }
  else
  {
    list_offset = 0;  // no mb aff or frame mb
    moving_block = Co_located->moving_block;
    co_located_mv = Co_located->mv;
    co_located_ref_idx = Co_located->ref_idx;
    co_located_ref_id = Co_located->ref_pic_id;
    max_y_cr = dec_picture->size_y_cr-1;
  }



  if (!img->MbaffFrameFlag)
  {
    for (l=0+list_offset;l<(2+list_offset);l++)
    {
      for(k = 0; k < listXsize[l]; k++)
      {
        listX[l][k]->chroma_vector_adjustment= 0;
        if(img->structure == TOP_FIELD && img->structure != listX[l][k]->structure)
          listX[l][k]->chroma_vector_adjustment = -2;
        if(img->structure == BOTTOM_FIELD && img->structure != listX[l][k]->structure)
          listX[l][k]->chroma_vector_adjustment = 2;
      }
    }
  }
  else
  {
    if (curr_mb_field)
    {
      for (l=0+list_offset;l<(2+list_offset);l++)
      {
        for(k = 0; k < listXsize[l]; k++)
        {
          listX[l][k]->chroma_vector_adjustment= 0;
          if(img->current_mb_nr % 2 == 0 && listX[l][k]->structure == BOTTOM_FIELD)
            listX[l][k]->chroma_vector_adjustment = -2;
          if(img->current_mb_nr % 2 == 1 && listX[l][k]->structure == TOP_FIELD)
            listX[l][k]->chroma_vector_adjustment = 2;
        }
      }
    }
    else
    {
      for (l=0+list_offset;l<(2+list_offset);l++)
      {
        for(k = 0; k < listXsize[l]; k++)
        {
          listX[l][k]->chroma_vector_adjustment= 0;
        }
      }
    }
  }

  mv_mul=4;
  f1=8;
  f2=7;

  f3=f1*f1;
  f4=f3/2;

  // luma decoding **************************************************

  // get prediction for INTRA_MB_16x16
  if (IS_NEWINTRA (currMB))
  {
    intrapred_luma_16x16(img, currMB->i16mode);//判断宏块是否采用16x16帧内预测，如果是，则通过函数
	                                          //intrapred_luma_16x16()解码。然后对一个16x16宏块中的
	                                          //每个4x4块分别解码（每个块的解码顺序见毕厚杰书的图6.49（P118）），
	                                            //直到所有块解码完成为止。
  }

  if (img->type==B_SLICE && img->direct_type && (IS_DIRECT (currMB) || 
    (IS_P8x8(currMB) && !(currMB->b8mode[0] && currMB->b8mode[1] && currMB->b8mode[2] && currMB->b8mode[3]))))
  {
    int fw_rFrameL, fw_rFrameU, fw_rFrameUL, fw_rFrameUR;
    int bw_rFrameL, bw_rFrameU, bw_rFrameUL, bw_rFrameUR;    
    
    PixelPos mb_left, mb_up, mb_upleft, mb_upright;              
    
    getLuma4x4Neighbour(img->current_mb_nr,0,0,-1, 0,&mb_left);
    getLuma4x4Neighbour(img->current_mb_nr,0,0, 0,-1,&mb_up);
    getLuma4x4Neighbour(img->current_mb_nr,0,0,16, -1,&mb_upright);
    getLuma4x4Neighbour(img->current_mb_nr,0,0, -1,-1,&mb_upleft);

    if (!img->MbaffFrameFlag)
    {
      fw_rFrameL = mb_left.available ? dec_picture->ref_idx[LIST_0][mb_left.pos_x][mb_left.pos_y] : -1;
      fw_rFrameU = mb_up.available ? dec_picture->ref_idx[LIST_0][mb_up.pos_x][mb_up.pos_y] : -1;
      fw_rFrameUL = mb_upleft.available ? dec_picture->ref_idx[LIST_0][mb_upleft.pos_x][mb_upleft.pos_y] : -1;
      fw_rFrameUR = mb_upright.available ? dec_picture->ref_idx[LIST_0][mb_upright.pos_x][mb_upright.pos_y] : fw_rFrameUL;      
      
      bw_rFrameL = mb_left.available ? dec_picture->ref_idx[LIST_1][mb_left.pos_x][mb_left.pos_y] : -1;
      bw_rFrameU = mb_up.available ? dec_picture->ref_idx[LIST_1][mb_up.pos_x][mb_up.pos_y] : -1;
      bw_rFrameUL = mb_upleft.available ? dec_picture->ref_idx[LIST_1][mb_upleft.pos_x][mb_upleft.pos_y] : -1;
      bw_rFrameUR = mb_upright.available ? dec_picture->ref_idx[LIST_1][mb_upright.pos_x][mb_upright.pos_y] : bw_rFrameUL;      
    }
    else
    {
      if (img->mb_data[img->current_mb_nr].mb_field)
      {
        fw_rFrameL = mb_left.available ? 
          img->mb_data[mb_left.mb_addr].mb_field  || dec_picture->ref_idx[LIST_0][mb_left.pos_x][mb_left.pos_y] < 0? 
          dec_picture->ref_idx[LIST_0][mb_left.pos_x][mb_left.pos_y] : 
          dec_picture->ref_idx[LIST_0][mb_left.pos_x][mb_left.pos_y] * 2: -1;

        fw_rFrameU = mb_up.available ? 
          img->mb_data[mb_up.mb_addr].mb_field || dec_picture->ref_idx[LIST_0][mb_up.pos_x][mb_up.pos_y] < 0? 
          dec_picture->ref_idx[LIST_0][mb_up.pos_x][mb_up.pos_y] : 
        dec_picture->ref_idx[LIST_0][mb_up.pos_x][mb_up.pos_y] * 2: -1;

        fw_rFrameUL = mb_upleft.available ? 
          img->mb_data[mb_upleft.mb_addr].mb_field || dec_picture->ref_idx[LIST_0][mb_upleft.pos_x][mb_upleft.pos_y] < 0?         
          dec_picture->ref_idx[LIST_0][mb_upleft.pos_x][mb_upleft.pos_y] : 
        dec_picture->ref_idx[LIST_0][mb_upleft.pos_x][mb_upleft.pos_y] *2: -1;      

        fw_rFrameUR = mb_upright.available ? 
          img->mb_data[mb_upright.mb_addr].mb_field || dec_picture->ref_idx[LIST_0][mb_upright.pos_x][mb_upright.pos_y] < 0?
          dec_picture->ref_idx[LIST_0][mb_upright.pos_x][mb_upright.pos_y] : 
        dec_picture->ref_idx[LIST_0][mb_upright.pos_x][mb_upright.pos_y] * 2: fw_rFrameUL;      
        
        bw_rFrameL = mb_left.available ? 
          img->mb_data[mb_left.mb_addr].mb_field || dec_picture->ref_idx[LIST_1][mb_left.pos_x][mb_left.pos_y] < 0? 
          dec_picture->ref_idx[LIST_1][mb_left.pos_x][mb_left.pos_y] : 
        dec_picture->ref_idx[LIST_1][mb_left.pos_x][mb_left.pos_y] * 2: -1;

        bw_rFrameU = mb_up.available ? 
          img->mb_data[mb_up.mb_addr].mb_field || dec_picture->ref_idx[LIST_1][mb_up.pos_x][mb_up.pos_y] < 0? 
          dec_picture->ref_idx[LIST_1][mb_up.pos_x][mb_up.pos_y] : 
        dec_picture->ref_idx[LIST_1][mb_up.pos_x][mb_up.pos_y] * 2: -1;

        bw_rFrameUL = mb_upleft.available ? 
          img->mb_data[mb_upleft.mb_addr].mb_field || dec_picture->ref_idx[LIST_1][mb_upleft.pos_x][mb_upleft.pos_y] < 0?         
          dec_picture->ref_idx[LIST_1][mb_upleft.pos_x][mb_upleft.pos_y] : 
        dec_picture->ref_idx[LIST_1][mb_upleft.pos_x][mb_upleft.pos_y] *2: -1;      

        bw_rFrameUR = mb_upright.available ? 
          img->mb_data[mb_upright.mb_addr].mb_field || dec_picture->ref_idx[LIST_1][mb_upright.pos_x][mb_upright.pos_y] < 0?         
          dec_picture->ref_idx[LIST_1][mb_upright.pos_x][mb_upright.pos_y] : 
        dec_picture->ref_idx[LIST_1][mb_upright.pos_x][mb_upright.pos_y] * 2: bw_rFrameUL;              
      }
      else
      {
        fw_rFrameL = mb_left.available ? 
          img->mb_data[mb_left.mb_addr].mb_field || dec_picture->ref_idx[LIST_0][mb_left.pos_x][mb_left.pos_y]  < 0 ?
          dec_picture->ref_idx[LIST_0][mb_left.pos_x][mb_left.pos_y] >> 1 : 
        dec_picture->ref_idx[LIST_0][mb_left.pos_x][mb_left.pos_y]: -1;
        
        fw_rFrameU = mb_up.available ? 
          img->mb_data[mb_up.mb_addr].mb_field || dec_picture->ref_idx[LIST_0][mb_up.pos_x][mb_up.pos_y] < 0 ?
          dec_picture->ref_idx[LIST_0][mb_up.pos_x][mb_up.pos_y] >> 1 :  
        dec_picture->ref_idx[LIST_0][mb_up.pos_x][mb_up.pos_y] : -1;
        
        fw_rFrameUL = mb_upleft.available ? 
          img->mb_data[mb_upleft.mb_addr].mb_field || dec_picture->ref_idx[LIST_0][mb_upleft.pos_x][mb_upleft.pos_y] < 0 ?
          dec_picture->ref_idx[LIST_0][mb_upleft.pos_x][mb_upleft.pos_y]>> 1 : 
        dec_picture->ref_idx[LIST_0][mb_upleft.pos_x][mb_upleft.pos_y] : -1;      
        
        fw_rFrameUR = mb_upright.available ? 
          img->mb_data[mb_upright.mb_addr].mb_field || dec_picture->ref_idx[LIST_0][mb_upright.pos_x][mb_upright.pos_y] < 0 ? 
          dec_picture->ref_idx[LIST_0][mb_upright.pos_x][mb_upright.pos_y] >> 1 :  
        dec_picture->ref_idx[LIST_0][mb_upright.pos_x][mb_upright.pos_y] : fw_rFrameUL;      
        
        bw_rFrameL = mb_left.available ? 
          img->mb_data[mb_left.mb_addr].mb_field || dec_picture->ref_idx[LIST_1][mb_left.pos_x][mb_left.pos_y] < 0 ?
          dec_picture->ref_idx[LIST_1][mb_left.pos_x][mb_left.pos_y] >> 1 :  
        dec_picture->ref_idx[LIST_1][mb_left.pos_x][mb_left.pos_y] : -1;
        
        bw_rFrameU = mb_up.available ? 
          img->mb_data[mb_up.mb_addr].mb_field || dec_picture->ref_idx[LIST_1][mb_up.pos_x][mb_up.pos_y] < 0 ?
          dec_picture->ref_idx[LIST_1][mb_up.pos_x][mb_up.pos_y] >> 1 : 
        dec_picture->ref_idx[LIST_1][mb_up.pos_x][mb_up.pos_y] : -1;
        
        bw_rFrameUL = mb_upleft.available ? 
          img->mb_data[mb_upleft.mb_addr].mb_field || dec_picture->ref_idx[LIST_1][mb_upleft.pos_x][mb_upleft.pos_y] < 0 ?
          dec_picture->ref_idx[LIST_1][mb_upleft.pos_x][mb_upleft.pos_y] >> 1 : 
        dec_picture->ref_idx[LIST_1][mb_upleft.pos_x][mb_upleft.pos_y] : -1;      
        
        bw_rFrameUR = mb_upright.available ? 
          img->mb_data[mb_upright.mb_addr].mb_field || dec_picture->ref_idx[LIST_1][mb_upright.pos_x][mb_upright.pos_y] < 0 ?
          dec_picture->ref_idx[LIST_1][mb_upright.pos_x][mb_upright.pos_y] >> 1: 
        dec_picture->ref_idx[LIST_1][mb_upright.pos_x][mb_upright.pos_y] : bw_rFrameUL;      
      }
    }
    

    fw_rFrame = (fw_rFrameL >= 0 && fw_rFrameU >= 0) ? min(fw_rFrameL,fw_rFrameU): max(fw_rFrameL,fw_rFrameU);
    fw_rFrame = (fw_rFrame >= 0 && fw_rFrameUR >= 0) ? min(fw_rFrame,fw_rFrameUR): max(fw_rFrame,fw_rFrameUR);
    
    bw_rFrame = (bw_rFrameL >= 0 && bw_rFrameU >= 0) ? min(bw_rFrameL,bw_rFrameU): max(bw_rFrameL,bw_rFrameU);
    bw_rFrame = (bw_rFrame >= 0 && bw_rFrameUR >= 0) ? min(bw_rFrame,bw_rFrameUR): max(bw_rFrame,bw_rFrameUR);        
    
    if (fw_rFrame >=0)
      SetMotionVectorPredictor (img, pmvfw, pmvfw+1, fw_rFrame, LIST_0, dec_picture->ref_idx, dec_picture->mv, 0, 0, 16, 16);
    
    if (bw_rFrame >=0)
      SetMotionVectorPredictor (img, pmvbw, pmvbw+1, bw_rFrame, LIST_1, dec_picture->ref_idx, dec_picture->mv, 0, 0, 16, 16);

  }

  for (block8x8=0; block8x8<4; block8x8++)	//++ 四个8*8亮度块的循环
  {
    for (k = block8x8*4; k < block8x8*4+4; k ++)	//++ 每个8*8亮度块内的四个4*4亮度块的循环
    {
      i = (decode_block_scan[k] & 3);			//++ 4*4亮度块在本宏块中的横坐标（以4*4块为单位）
      j = ((decode_block_scan[k] >> 2) & 3);	//++ 4*4亮度块在本宏块中的纵坐标（以4*4块为单位）
      
      ioff=i*4;				//++ 4*4亮度块左上角像素在本宏块中的横坐标
      i4=img->block_x+i;	//++ 4*4亮度块在本帧中的横坐标（以4*4块为单位）
      
      joff=j*4;				//++ 4*4亮度块左上角像素在本宏块中的纵坐标
      j4=img->block_y+j;	//++ 4*4亮度块在本帧中的纵坐标（以4*4块为单位）

      mv_mode  = currMB->b8mode[2*(j/2)+(i/2)];
      pred_dir = currMB->b8pdir[2*(j/2)+(i/2)];
      
      assert (pred_dir<=2);

      // PREDICTION
      if (mv_mode==IBLOCK)
      {
        //===== INTRA PREDICTION =====
		//++ 4*4帧内预测
        if (intrapred(img,ioff,joff,i4,j4)==SEARCH_SYNC)  /* make 4x4 prediction block mpr from given prediction img->mb_mode */
          return SEARCH_SYNC;                   /* bit error */
      }
      else if (!IS_NEWINTRA (currMB))
      {
        if (pred_dir != 2)
        {
          //===== FORWARD/BACKWARD PREDICTION =====
          fw_refframe = ref_idx  = dec_picture->ref_idx[LIST_0 + pred_dir][i4][j4];
          mv_array = dec_picture->mv[LIST_0 + pred_dir];
          list     = listX[0+list_offset+ pred_dir];
          vec1_x = i4*4*mv_mul + mv_array[i4][j4][0];	//++ 当前4*4块的参考块左上角像素在本图像帧中的坐标（以1/4像素距离为单位）

          if (!curr_mb_field)
          {
            vec1_y = j4*4*mv_mul + mv_array[i4][j4][1];	//++ 当前4*4块的参考块左上角像素在本图像帧中的坐标（以1/4像素距离为单位）
          }
          else
          {
            if (mb_nr%2 == 0) 
              vec1_y = (img->block_y * 2 + joff) * mv_mul + mv_array[i4][j4][1];	//++ 当前4*4块的参考块左上角像素在本图像帧中的坐标（以1/4像素距离为单位）
            else
              vec1_y = ((img->block_y-4) * 2 + joff)* mv_mul + mv_array[i4][j4][1];	//++ 当前4*4块的参考块左上角像素在本图像帧中的坐标（以1/4像素距离为单位）
          }

          get_block (ref_idx, list, vec1_x, vec1_y, img, tmp_block);

          if (img->apply_weights)
          {
            if (((active_pps->weighted_pred_flag&&(img->type==P_SLICE|| img->type == SP_SLICE))||
                (active_pps->weighted_bipred_idc==1 && (img->type==B_SLICE))) && curr_mb_field)
            {
              ref_idx >>=1;
            }

            for(ii=0;ii<BLOCK_SIZE;ii++)
              for(jj=0;jj<BLOCK_SIZE;jj++)  
                img->mpr[ii+ioff][jj+joff] = Clip1(((img->wp_weight[pred_dir][ref_idx][0] *  tmp_block[ii][jj]+ img->wp_round_luma) >>img->luma_log2_weight_denom)  + img->wp_offset[pred_dir][fw_refframe>>curr_mb_field][0] );
          }
          else
          {
            for(ii=0;ii<BLOCK_SIZE;ii++)
              for(jj=0;jj<BLOCK_SIZE;jj++)  
                img->mpr[ii+ioff][jj+joff] = tmp_block[ii][jj];
          }
        }
        else
        {
          if (mv_mode != 0)
          {
            //===== BI-DIRECTIONAL PREDICTION =====
            fw_mv_array = dec_picture->mv[LIST_0];
            bw_mv_array = dec_picture->mv[LIST_1];
            
            fw_refframe = dec_picture->ref_idx[LIST_0][i4][j4];
            bw_refframe = dec_picture->ref_idx[LIST_1][i4][j4];
            fw_ref_idx = fw_refframe;
            bw_ref_idx = bw_refframe;
          }
          else
          {
            //===== DIRECT PREDICTION =====
            fw_mv_array = dec_picture->mv[LIST_0];	//++ 前向运动矢量数组
            bw_mv_array = dec_picture->mv[LIST_1];	//++ 后向运动矢量数组
            bw_refframe = 0;

            if (img->direct_type )	//++ 空间直接模式
            {
              int imgblock_y= ((img->MbaffFrameFlag)&&(currMB->mb_field))? (img->current_mb_nr%2) ? (img->block_y-4)/2:img->block_y/2: img->block_y;
              int j6=imgblock_y+j;


              if (fw_rFrame >=0)
              {
                if (!fw_rFrame  && ((!moving_block[i4][j6]) && (!listX[1+list_offset][0]->is_long_term)))
                {
                  dec_picture->mv  [LIST_0][i4][j4][0]= 0;
                  dec_picture->mv  [LIST_0][i4][j4][1]= 0;
                  dec_picture->ref_idx[LIST_0][i4][j4] = 0;
                }
                else
                {
                  dec_picture->mv  [LIST_0][i4][j4][0]= pmvfw[0];
                  dec_picture->mv  [LIST_0][i4][j4][1]= pmvfw[1];
                  dec_picture->ref_idx[LIST_0][i4][j4] = fw_rFrame;
                }
              }
              else
              {
                dec_picture->ref_idx[LIST_0][i4][j4] = -1;                
                dec_picture->mv  [LIST_0][i4][j4][0]= 0;
                dec_picture->mv  [LIST_0][i4][j4][1]= 0;
              }
              
              if (bw_rFrame >=0)
              {
                if  (bw_rFrame==0 && ((!moving_block[i4][j6]) && (!listX[1+list_offset][0]->is_long_term)))
                {                  
                  
                  dec_picture->mv  [LIST_1][i4][j4][0]= 0;
                  dec_picture->mv  [LIST_1][i4][j4][1]= 0;
                  dec_picture->ref_idx[LIST_1][i4][j4] = bw_rFrame;
                  
                }
                else
                {
                  dec_picture->mv  [LIST_1][i4][j4][0]= pmvbw[0];
                  dec_picture->mv  [LIST_1][i4][j4][1]= pmvbw[1];
                  
                  dec_picture->ref_idx[LIST_1][i4][j4] = bw_rFrame;
                }               
              }
              else
              {                  
                dec_picture->mv  [LIST_1][i4][j4][0]=0;
                dec_picture->mv  [LIST_1][i4][j4][1]=0;
                dec_picture->ref_idx[LIST_1][i4][j4] = -1;
                
              }
              
              if (fw_rFrame < 0 && bw_rFrame < 0)
              {
                dec_picture->ref_idx[LIST_0][i4][j4] = 0;
                dec_picture->ref_idx[LIST_1][i4][j4] = 0;
              }
              
              fw_refframe = (dec_picture->ref_idx[LIST_0][i4][j4]!=-1) ? dec_picture->ref_idx[LIST_0][i4][j4]:0;
              bw_refframe = (dec_picture->ref_idx[LIST_1][i4][j4]!=-1) ? dec_picture->ref_idx[LIST_1][i4][j4]:0;
              
              fw_ref_idx = fw_refframe;
              bw_ref_idx = bw_refframe;
              
              if      (dec_picture->ref_idx[LIST_1][i4][j4]==-1) direct_pdir = 0;
              else if (dec_picture->ref_idx[LIST_0][i4][j4]==-1) direct_pdir = 1;
              else                                               direct_pdir = 2;
              
            }
            else // Temporal Mode	//++ 时间直接模式
            {
              
              int imgblock_y= ((img->MbaffFrameFlag)&&(currMB->mb_field))? (img->current_mb_nr%2) ? (img->block_y-4)/2:img->block_y/2: img->block_y;
              int j6= imgblock_y + j;
              
              int refList = (co_located_ref_idx[LIST_0][i4][j6]== -1 ? LIST_1 : LIST_0);	//++ 如果“相对位置块”的前向参考帧不可用则采用其后向参考帧
              int ref_idx =  co_located_ref_idx[refList][i4][j6];


              if(ref_idx==-1) // co-located is intra mode
              {
                for(hv=0; hv<2; hv++)
                {
                  dec_picture->mv  [LIST_0][i4][j4][hv]=0;
                  dec_picture->mv  [LIST_1][i4][j4][hv]=0;                    
                }
                
                dec_picture->ref_idx[LIST_0][i4][j4] = 0;
                dec_picture->ref_idx[LIST_1][i4][j4] = 0;
                
                fw_refframe = 0;
                fw_ref_idx = 0;
              }
              else // co-located skip or inter mode
              {
                int mapped_idx=0;
                int iref;          
                
                {
                  for (iref=0;iref<min(img->num_ref_idx_l0_active,listXsize[LIST_0 + list_offset]);iref++)
                  {
                    if (dec_picture->ref_pic_num[img->current_slice_nr][LIST_0 + list_offset][iref]==co_located_ref_id[refList][i4][j6])
                    {
                      mapped_idx=iref;
                      break;
                    }
                    else //! invalid index. Default to zero even though this case should not happen
                    {                        
                      mapped_idx=INVALIDINDEX;
                    }
                  }
                  if (INVALIDINDEX == mapped_idx)
                  {   
                    error("temporal direct error\ncolocated block has ref that is unavailable",-1111);
                  }
                }
          
                fw_ref_idx = mapped_idx;
                mv_scale = img->mvscale[LIST_0 + list_offset][mapped_idx];

                //! In such case, an array is needed for each different reference.
//                if (mv_scale == 9999 || Co_located->is_long_term)
                if (mv_scale == 9999 || listX[LIST_0+list_offset][mapped_idx]->is_long_term)
                {
                  dec_picture->mv  [LIST_0][i4][j4][0]=co_located_mv[refList][i4][j6][0];
                  dec_picture->mv  [LIST_0][i4][j4][1]=co_located_mv[refList][i4][j6][1];

                  dec_picture->mv  [LIST_1][i4][j4][0]=0;
                  dec_picture->mv  [LIST_1][i4][j4][1]=0;
                }
                else
                {
                  dec_picture->mv  [LIST_0][i4][j4][0]=(mv_scale * co_located_mv[refList][i4][j6][0] + 128 ) >> 8;	//++ 参见标准文档公式8-196
                  dec_picture->mv  [LIST_0][i4][j4][1]=(mv_scale * co_located_mv[refList][i4][j6][1] + 128 ) >> 8;	//++ 参见标准文档公式8-196
                  
                  dec_picture->mv  [LIST_1][i4][j4][0]=dec_picture->mv  [LIST_0][i4][j4][0] - co_located_mv[refList][i4][j6][0] ;
                  dec_picture->mv  [LIST_1][i4][j4][1]=dec_picture->mv  [LIST_0][i4][j4][1] - co_located_mv[refList][i4][j6][1] ;
                }
                
                fw_refframe = dec_picture->ref_idx[LIST_0][i4][j4] = mapped_idx; //listX[1][0]->ref_idx[refList][i4][j4];
                bw_refframe = dec_picture->ref_idx[LIST_1][i4][j4] = 0;
                
                fw_ref_idx = fw_refframe;
                bw_ref_idx = bw_refframe;
              }
            }
            // store reference picture ID determined by direct mode
            dec_picture->ref_pic_id[LIST_0][i4][j4] = dec_picture->ref_pic_num[img->current_slice_nr][LIST_0 + list_offset][dec_picture->ref_idx[LIST_0][i4][j4]];
            dec_picture->ref_pic_id[LIST_1][i4][j4] = dec_picture->ref_pic_num[img->current_slice_nr][LIST_1 + list_offset][dec_picture->ref_idx[LIST_1][i4][j4]];  
          }
                 
          if (mv_mode==0 && img->direct_type )
          {
            if (dec_picture->ref_idx[LIST_0][i4][j4] >= 0)
            {
              
              vec1_x = i4*4*mv_mul + fw_mv_array[i4][j4][0];
              if (!curr_mb_field)
              {
                vec1_y = j4*4*mv_mul + fw_mv_array[i4][j4][1];
              }
              else
              {
                if (mb_nr%2 == 0)
                {
                  vec1_y = (img->block_y * 2 + joff) * mv_mul + fw_mv_array[i4][j4][1];
                }
                else
                {
                  vec1_y = ((img->block_y-4) * 2 + joff)* mv_mul + fw_mv_array[i4][j4][1];
                }
              }               
              get_block(fw_refframe, listX[0+list_offset], vec1_x, vec1_y, img, tmp_block);
            }
                  
            if (dec_picture->ref_idx[LIST_1][i4][j4] >= 0)
            {
              vec2_x = i4*4*mv_mul + bw_mv_array[i4][j4][0];
              if (!curr_mb_field)
              {
                vec2_y = j4*4*mv_mul + bw_mv_array[i4][j4][1];
              }
              else
              {
                if (mb_nr%2 == 0)
                {
                  vec2_y = (img->block_y * 2 + joff) * mv_mul + bw_mv_array[i4][j4][1];
                }
                else
                {
                  vec2_y = ((img->block_y-4) * 2 + joff)* mv_mul + bw_mv_array[i4][j4][1];
                }            
              }
              get_block(bw_refframe, listX[1+list_offset], vec2_x, vec2_y, img, tmp_blockbw);
            }              
          }
          else
          {
            vec1_x = i4*4*mv_mul + fw_mv_array[i4][j4][0];
            vec2_x = i4*4*mv_mul + bw_mv_array[i4][j4][0];
            
            if (!curr_mb_field)
            {
              vec1_y = j4*4*mv_mul + fw_mv_array[i4][j4][1];
              vec2_y = j4*4*mv_mul + bw_mv_array[i4][j4][1];
            }
            else
            {
              if (mb_nr%2 == 0)
              {
                vec1_y = (img->block_y * 2 + joff) * mv_mul + fw_mv_array[i4][j4][1];
                vec2_y = (img->block_y * 2 + joff) * mv_mul + bw_mv_array[i4][j4][1];
              }
              else
              {
                vec1_y = ((img->block_y-4) * 2 + joff)* mv_mul + fw_mv_array[i4][j4][1];
                vec2_y = ((img->block_y-4) * 2 + joff)* mv_mul + bw_mv_array[i4][j4][1];
              }
            }
            
            get_block(fw_refframe, listX[0+list_offset], vec1_x, vec1_y, img, tmp_block);
            get_block(bw_refframe, listX[1+list_offset], vec2_x, vec2_y, img, tmp_blockbw);
          }
          
          if (mv_mode==0 && img->direct_type && direct_pdir==0)
          {
            if (img->apply_weights)
            {
              if (((active_pps->weighted_pred_flag&&(img->type==P_SLICE|| img->type == SP_SLICE))||
                (active_pps->weighted_bipred_idc==1 && (img->type==B_SLICE))) && curr_mb_field)
              {
                fw_ref_idx >>=1;
              }
              for(ii=0;ii<BLOCK_SIZE;ii++)
                for(jj=0;jj<BLOCK_SIZE;jj++)  
                  img->mpr[ii+ioff][jj+joff] = Clip1(((tmp_block[ii][jj] * img->wp_weight[0][fw_ref_idx][0]  
                                               + img->wp_round_luma)>>img->luma_log2_weight_denom) 
                                               + img->wp_offset[0][fw_refframe>>curr_mb_field][0]);
            }
            else
            {
              for(ii=0;ii<BLOCK_SIZE;ii++)
                for(jj=0;jj<BLOCK_SIZE;jj++)  
                  img->mpr[ii+ioff][jj+joff] = tmp_block[ii][jj];
            }
          }
          else if (mv_mode==0 && img->direct_type && direct_pdir==1)
          {              
            if (img->apply_weights)
            {
              if (((active_pps->weighted_pred_flag&&(img->type==P_SLICE|| img->type == SP_SLICE))||
                (active_pps->weighted_bipred_idc==1 && (img->type==B_SLICE))) && curr_mb_field)
              {
                fw_ref_idx >>=1;
                bw_ref_idx >>=1;
              }

              for(ii=0;ii<BLOCK_SIZE;ii++)
                for(jj=0;jj<BLOCK_SIZE;jj++)  
                  img->mpr[ii+ioff][jj+joff] = Clip1(((tmp_blockbw[ii][jj] * img->wp_weight[1][bw_ref_idx][0] 
                                               + img->wp_round_luma)>>img->luma_log2_weight_denom) 
                                               + img->wp_offset[1][bw_refframe>>curr_mb_field][0]);
            }
            else
            {
              for(ii=0;ii<BLOCK_SIZE;ii++)
                for(jj=0;jj<BLOCK_SIZE;jj++)  
                  img->mpr[ii+ioff][jj+joff] = tmp_blockbw[ii][jj];
            }
          }
          else if(img->apply_weights)
          {
            int alpha_fw, alpha_bw;
            int wt_list_offset = (active_pps->weighted_bipred_idc==2)?list_offset:0;

            if (mv_mode==0 && img->direct_type==0 )bw_ref_idx=0;    //temporal direct 

            if (((active_pps->weighted_pred_flag&&(img->type==P_SLICE|| img->type == SP_SLICE))||
              (active_pps->weighted_bipred_idc==1 && (img->type==B_SLICE))) && curr_mb_field)
            {
              fw_ref_idx >>=1;
              bw_ref_idx >>=1;
            }

            alpha_fw = img->wbp_weight[0+wt_list_offset][fw_ref_idx][bw_ref_idx][0];
            alpha_bw = img->wbp_weight[1+wt_list_offset][fw_ref_idx][bw_ref_idx][0];

            for(ii=0;ii<BLOCK_SIZE;ii++)
              for(jj=0;jj<BLOCK_SIZE;jj++)  
                img->mpr[ii+ioff][jj+joff] = (int)Clip1(((alpha_fw * tmp_block[ii][jj] + alpha_bw * tmp_blockbw[ii][jj]  
                                             + (1<<img->luma_log2_weight_denom)) >> (img->luma_log2_weight_denom+1)) 
                                             + ((img->wp_offset[wt_list_offset+0][fw_ref_idx][0] + img->wp_offset[wt_list_offset+1][bw_ref_idx][0] + 1) >>1));
          }
          else
          {
            for(ii=0;ii<BLOCK_SIZE;ii++)
              for(jj=0;jj<BLOCK_SIZE;jj++)  
                img->mpr[ii+ioff][jj+joff] = (tmp_block[ii][jj]+tmp_blockbw[ii][jj]+1)/2;
          }
        }
      }
      if (smb && mv_mode!=IBLOCK)
      {
        itrans_sp(img,ioff,joff,i,j);
      }
      else
      {	//++ 反变换
        itrans   (img,ioff,joff,i,j);      // use DCT transform and make 4x4 block m7 from prediction block mpr
      }

      for(ii=0;ii<BLOCK_SIZE;ii++)
      {     
        for(jj=0;jj<BLOCK_SIZE;jj++)
        {
          dec_picture->imgY[j4*BLOCK_SIZE+jj][i4*BLOCK_SIZE+ii]=img->m7[ii][jj]; // contruct picture from 4x4 blocks
        }
      }
    }
  }
  
  // chroma decoding *******************************************************
  for(uv=0;uv<2;uv++)
  {
    intra_prediction = IS_INTRA (currMB);
    
    if (intra_prediction)
    {
      intrapred_chroma(img, uv);
    }
    
    for (j=4;j<6;j++)			//++ 四个4*4色度块循环的第一层循环
    {
      joff=(j-4)*4;				//++ 当前4*4块左上角像素在在本色度块内的纵坐标
      j4=img->pix_c_y+joff;		//++ 当前4*4块左上角像素在本帧内的纵坐标
      for(i=0;i<2;i++)			//++ 四个4*4色度块的循环的第二层循环
      {
        ioff=i*4;				//++ 当前4*4块左上角像素在在本色度块内的横坐标
        i4=img->pix_c_x+ioff;	//++ 当前4*4块左上角像素在本帧内的横坐标
        
        mv_mode  = currMB->b8mode[2*(j-4)+i];
        pred_dir = currMB->b8pdir[2*(j-4)+i];
        assert (pred_dir<=2);

        if (!intra_prediction)
        {
          if (pred_dir != 2)
          {
            //--- FORWARD/BACKWARD PREDICTION ---
            mv_array = dec_picture->mv[LIST_0 + pred_dir];
            list = listX[0+list_offset+pred_dir];
            for(jj=0;jj<4;jj++)		//++ 每个4*4块内的16个像素的第一层循环
            {
              jf=(j4+jj)/2;			//++ 当前像素在本帧内的纵坐标（以2*2块为单位），主要是为了计算当前像素对应哪一个4*4的亮度块（因为当前像素与对应的4*4亮度块采用同样的参考帧）
              for(ii=0;ii<4;ii++)	//++ 每个4*4块内的16个像素的第二层循环
              {
                if1=(i4+ii)/2;		//++ 当前像素在本帧内的横坐标（以2*2块为单位），主要是为了计算当前像素对应哪一个4*4的亮度块（因为当前像素与对应的4*4亮度块采用同样的参考帧）
                
                fw_refframe = ref_idx   = dec_picture->ref_idx[LIST_0+pred_dir][if1][jf];	//++ 当前像素的参考帧
                
                i1=(img->pix_c_x+ii+ioff)*f1+mv_array[if1][jf][0];	//++ 当前像素的运动向量X分量（1/4像素为单位）

                if (!curr_mb_field)
                  j1=(img->pix_c_y+jj+joff)*f1+mv_array[if1][jf][1];	//++ 当前像素的运动向量Y分量（1/4像素为单位）
                else
                {
                  if (mb_nr%2 == 0) 
                    j1=(img->pix_c_y)/2*f1 + (jj+joff)*f1+mv_array[if1][jf][1];
                  else
                    j1=(img->pix_c_y-8)/2*f1 + (jj+joff)*f1 +mv_array[if1][jf][1];
                }
                
                j1 += list[ref_idx]->chroma_vector_adjustment;
                
                ii0=max (0, min (i1>>3, img->width_cr-1));
                jj0=max (0, min (j1>>3, max_y_cr));
                ii1=max (0, min ((i1>>3)+1, img->width_cr-1));
                jj1=max (0, min ((j1>>3)+1, max_y_cr));
                
                if1=(i1 & f2);	//++ 当前像素的参考像素（不一定在整像素位置）与其左方最邻近整像素的距离（以1/4像素距离为单位）
                jf1=(j1 & f2);	//++ 当前像素的参考像素（不一定在整像素位置）与其上方最邻近整像素的距离（以1/4像素距离为单位）
                if0=f1-if1;		//++ 当前像素的参考像素（不一定在整像素位置）与其右方最邻近整像素的距离（以1/4像素距离为单位）
                jf0=f1-jf1;		//++ 当前像素的参考像素（不一定在整像素位置）与其下方最邻近整像素的距离（以1/4像素距离为单位）
                
                if (img->apply_weights)
                {
                  pred = (if0*jf0*list[ref_idx]->imgUV[uv][jj0][ii0]+
                          if1*jf0*list[ref_idx]->imgUV[uv][jj0][ii1]+
                          if0*jf1*list[ref_idx]->imgUV[uv][jj1][ii0]+
                          if1*jf1*list[ref_idx]->imgUV[uv][jj1][ii1]+f4)>>6;
                  if (((active_pps->weighted_pred_flag&&(img->type==P_SLICE|| img->type == SP_SLICE))||
                    (active_pps->weighted_bipred_idc==1 && (img->type==B_SLICE))) && curr_mb_field)
                  {
                    ref_idx >>=1;
                  }

                  img->mpr[ii+ioff][jj+joff] = Clip1(((img->wp_weight[pred_dir][ref_idx][uv+1] * pred  + img->wp_round_chroma)>>img->chroma_log2_weight_denom) + img->wp_offset[pred_dir][ref_idx][uv+1]);
                }
                else
                {
                  img->mpr[ii+ioff][jj+joff]=(if0*jf0*list[ref_idx]->imgUV[uv][jj0][ii0]+
                                              if1*jf0*list[ref_idx]->imgUV[uv][jj0][ii1]+
                                              if0*jf1*list[ref_idx]->imgUV[uv][jj1][ii0]+
                                              if1*jf1*list[ref_idx]->imgUV[uv][jj1][ii1]+f4)>>6;
                }
              }
            }
          }
          else
          {
            if (mv_mode != 0)
            {
              //===== BI-DIRECTIONAL PREDICTION =====
              fw_mv_array = dec_picture->mv[LIST_0];
              bw_mv_array = dec_picture->mv[LIST_1];
            }
            else
            {
              //===== DIRECT PREDICTION =====
              fw_mv_array = dec_picture->mv[LIST_0];
              bw_mv_array = dec_picture->mv[LIST_1];
            }
            
            for(jj=0;jj<4;jj++)
            {
              jf=(j4+jj)/2;
              for(ii=0;ii<4;ii++)
              {
                ifx=(i4+ii)/2;
                
                direct_pdir = 2;
                
                if (mv_mode != 0)
                {
                  fw_refframe = dec_picture->ref_idx[LIST_0][ifx][jf];
                  bw_refframe = dec_picture->ref_idx[LIST_1][ifx][jf];
                  
                  fw_ref_idx = fw_refframe;
                  bw_ref_idx = bw_refframe;
                }
                else
                {
                  if (!mv_mode && img->direct_type )	//++ 色度空间直接模式
                  {
                    if (dec_picture->ref_idx[LIST_0][(ifx/2)*2][2*(jf/2)]!=-1)
                    {
                      fw_refframe = dec_picture->ref_idx[LIST_0][(ifx/2)*2][2*(jf/2)];
                      fw_ref_idx = fw_refframe;
                    }
                    if (dec_picture->ref_idx[LIST_1][(ifx/2)*2][2*(jf/2)]!=-1)
                    {
                      bw_refframe = dec_picture->ref_idx[LIST_1][(ifx/2)*2][2*(jf/2)];
                      bw_ref_idx = bw_refframe;
                    }
                    
                    if      (dec_picture->ref_idx[LIST_1][(ifx/2)*2][2*(jf/2)]==-1) direct_pdir = 0;
                    else if (dec_picture->ref_idx[LIST_0][(ifx/2)*2][2*(jf/2)]==-1) direct_pdir = 1;
                    else                                                            direct_pdir = 2;
                  }                
                  else	//++ 色度时间直接模式
                  {
                    
                    fw_refframe = dec_picture->ref_idx[LIST_0][ifx][jf];
                    bw_refframe = dec_picture->ref_idx[LIST_1][ifx][jf];
                    
                    fw_ref_idx = fw_refframe;
                    bw_ref_idx = bw_refframe;
                  }
                }
                
                
                if (mv_mode==0 && img->direct_type )	//++ 色度空间直接模式
                {
                  if (direct_pdir == 0 || direct_pdir == 2)
                  {
                    i1=(img->pix_c_x+ii+ioff)*f1+fw_mv_array[ifx][jf][0];                  
                    
                    if (!curr_mb_field)
                    {
                      j1=(img->pix_c_y+jj+joff)*f1+fw_mv_array[ifx][jf][1];
                    }
                    else
                    {
                      if (mb_nr%2 == 0) 
                        j1=(img->pix_c_y)/2*f1 + (jj+joff)*f1+fw_mv_array[ifx][jf][1];
                      else
                        j1=(img->pix_c_y-8)/2*f1 + (jj+joff)*f1 +fw_mv_array[ifx][jf][1];
                    }
                
                    j1 += listX[0+list_offset][fw_refframe]->chroma_vector_adjustment;
                    
                    ii0=max (0, min (i1>>3, img->width_cr-1));
                    jj0=max (0, min (j1>>3, max_y_cr));
                    ii1=max (0, min ((i1>>3)+1, img->width_cr-1));
                    jj1=max (0, min ((j1>>3)+1, max_y_cr));
                    
                    if1=(i1 & f2);
                    jf1=(j1 & f2);
                    if0=f1-if1;
                    jf0=f1-jf1;
                    
                    fw_pred=(if0*jf0*listX[0+list_offset][fw_refframe]->imgUV[uv][jj0][ii0]+
                             if1*jf0*listX[0+list_offset][fw_refframe]->imgUV[uv][jj0][ii1]+
                             if0*jf1*listX[0+list_offset][fw_refframe]->imgUV[uv][jj1][ii0]+
                             if1*jf1*listX[0+list_offset][fw_refframe]->imgUV[uv][jj1][ii1]+f4)>>6;
                  }
                  if (direct_pdir == 1 || direct_pdir == 2)
                  {
                    i1=(img->pix_c_x+ii+ioff)*f1+bw_mv_array[ifx][jf][0];
                    
                    if (!curr_mb_field)
                    {
                      j1=(img->pix_c_y+jj+joff)*f1+bw_mv_array[ifx][jf][1];
                    }
                    else
                    {
                      if (mb_nr%2 == 0) 
                        j1=(img->pix_c_y)/2*f1 + (jj+joff)*f1+bw_mv_array[ifx][jf][1];
                      else
                        j1=(img->pix_c_y-8)/2*f1 + (jj+joff)*f1 +bw_mv_array[ifx][jf][1];
                    }
                    j1 += listX[1+list_offset][bw_refframe]->chroma_vector_adjustment;
                    
                    ii0=max (0, min (i1>>3, img->width_cr-1));
                    jj0=max (0, min (j1>>3, max_y_cr));
                    ii1=max (0, min ((i1>>3)+1, img->width_cr-1));
                    jj1=max (0, min ((j1>>3)+1, max_y_cr));
                    
                    if1=(i1 & f2);
                    jf1=(j1 & f2);
                    if0=f1-if1;
                    jf0=f1-jf1;
                
                    bw_pred=(if0*jf0*listX[1+list_offset][bw_refframe]->imgUV[uv][jj0][ii0]+
                             if1*jf0*listX[1+list_offset][bw_refframe]->imgUV[uv][jj0][ii1]+
                             if0*jf1*listX[1+list_offset][bw_refframe]->imgUV[uv][jj1][ii0]+
                             if1*jf1*listX[1+list_offset][bw_refframe]->imgUV[uv][jj1][ii1]+f4)>>6;
                  }

                }
                else	//++ 色度时间直接模式
                {
                  i1=(img->pix_c_x+ii+ioff)*f1+fw_mv_array[ifx][jf][0];

                  if (!curr_mb_field)
                  {
                    j1=(img->pix_c_y+jj+joff)*f1+fw_mv_array[ifx][jf][1];
                  }
                  else
                  {
                    if (mb_nr%2 == 0) 
                      j1=(img->pix_c_y)/2*f1 + (jj+joff)*f1+fw_mv_array[ifx][jf][1];
                    else
                      j1=(img->pix_c_y-8)/2*f1 + (jj+joff)*f1 +fw_mv_array[ifx][jf][1];
                  }

                  j1 += listX[0+list_offset][fw_refframe]->chroma_vector_adjustment;
                  
                  ii0=max (0, min (i1>>3, img->width_cr-1));
                  jj0=max (0, min (j1>>3, max_y_cr));
                  ii1=max (0, min ((i1>>3)+1, img->width_cr-1));
                  jj1=max (0, min ((j1>>3)+1, max_y_cr));
                  
                  if1=(i1 & f2);
                  jf1=(j1 & f2);
                  if0=f1-if1;
                  jf0=f1-jf1;
                  
                  fw_pred=(if0*jf0*listX[0+list_offset][fw_refframe]->imgUV[uv][jj0][ii0]+
                           if1*jf0*listX[0+list_offset][fw_refframe]->imgUV[uv][jj0][ii1]+
                           if0*jf1*listX[0+list_offset][fw_refframe]->imgUV[uv][jj1][ii0]+
                           if1*jf1*listX[0+list_offset][fw_refframe]->imgUV[uv][jj1][ii1]+f4)>>6;
                  
                  i1=(img->pix_c_x+ii+ioff)*f1+bw_mv_array[ifx][jf][0];
                  
                  if (!curr_mb_field)
                  {
                    j1=(img->pix_c_y+jj+joff)*f1+bw_mv_array[ifx][jf][1];
                  }
                  else
                  {
                    if (mb_nr%2 == 0) 
                      j1=(img->pix_c_y)/2*f1 + (jj+joff)*f1+bw_mv_array[ifx][jf][1];
                    else
                      j1=(img->pix_c_y-8)/2*f1 + (jj+joff)*f1 +bw_mv_array[ifx][jf][1];
                  }

                  j1 += listX[1+list_offset][bw_refframe]->chroma_vector_adjustment;

                  ii0=max (0, min (i1>>3, img->width_cr-1));
                  jj0=max (0, min (j1>>3, max_y_cr));
                  ii1=max (0, min ((i1>>3)+1, img->width_cr-1));
                  jj1=max (0, min ((j1>>3)+1, max_y_cr));
                  
                  if1=(i1 & f2);
                  jf1=(j1 & f2);
                  if0=f1-if1;
                  jf0=f1-jf1;
                  
                  bw_pred=(if0*jf0*listX[1+list_offset][bw_refframe]->imgUV[uv][jj0][ii0]+
                           if1*jf0*listX[1+list_offset][bw_refframe]->imgUV[uv][jj0][ii1]+
                           if0*jf1*listX[1+list_offset][bw_refframe]->imgUV[uv][jj1][ii0]+
                           if1*jf1*listX[1+list_offset][bw_refframe]->imgUV[uv][jj1][ii1]+f4)>>6;
                  
                }
                
                if (img->apply_weights)
                {
                  if (((active_pps->weighted_pred_flag&&(img->type==P_SLICE|| img->type == SP_SLICE))||
                    (active_pps->weighted_bipred_idc==1 && (img->type==B_SLICE))) && curr_mb_field)
                  {
                    fw_ref_idx >>=1;
                    bw_ref_idx >>=1;
                  }

                  if (img->direct_type && direct_pdir==1)
                  {
                    img->mpr[ii+ioff][jj+joff]= Clip1(((img->wp_weight[1][bw_ref_idx][uv+1] * bw_pred  + img->wp_round_chroma)>>img->chroma_log2_weight_denom) + img->wp_offset[1][bw_refframe>>curr_mb_field][uv+1]);   // Replaced with integer only operations
                  }
                  else if (img->direct_type && direct_pdir==0)
                  {
                    img->mpr[ii+ioff][jj+joff]=Clip1(((img->wp_weight[0][fw_ref_idx][uv+1] * fw_pred + img->wp_round_chroma)>>img->chroma_log2_weight_denom) + img->wp_offset[0][fw_refframe>>curr_mb_field][uv+1]);   // Replaced with integer only operations
                  }
                  else
                  {
                    
                    int wt_list_offset = (active_pps->weighted_bipred_idc==2)?list_offset:0;

                    int alpha_fw = img->wbp_weight[0+wt_list_offset][fw_ref_idx][bw_ref_idx][uv+1];
                    int alpha_bw = img->wbp_weight[1+wt_list_offset][fw_ref_idx][bw_ref_idx][uv+1];

                    img->mpr[ii+ioff][jj+joff]= Clip1(((alpha_fw * fw_pred + alpha_bw * bw_pred  + (1<<img->chroma_log2_weight_denom)) >> (img->chroma_log2_weight_denom + 1))+ ((img->wp_offset[wt_list_offset + 0][fw_ref_idx][uv+1] + img->wp_offset[wt_list_offset + 1][bw_ref_idx][uv+1] + 1)>>1) );
                  }
                }
                else
                {
                  if (img->direct_type && direct_pdir==1)
                  {
                    img->mpr[ii+ioff][jj+joff]=bw_pred;
                  }
                  else if (img->direct_type && direct_pdir==0)
                  {
                    img->mpr[ii+ioff][jj+joff]=fw_pred;
                  } 
                  else
                  {
                    img->mpr[ii+ioff][jj+joff]=(fw_pred + bw_pred + 1 )/2;
                  }
                }
              }
            }
          }
        }
        
        if (!smb)
        {
          itrans(img,ioff,joff,2*uv+i,j);
          for(ii=0;ii<4;ii++)
            for(jj=0;jj<4;jj++)
            {
              dec_picture->imgUV[uv][j4+jj][i4+ii]=img->m7[ii][jj];
            }
        }
      }
    }
    
    if(smb)
    {
      itrans_sp_chroma(img,2*uv);
      for (j=4;j<6;j++)
      {
        joff=(j-4)*4;
        j4=img->pix_c_y+joff;
        for(i=0;i<2;i++)
        {
          ioff=i*4;
          i4=img->pix_c_x+ioff;
          itrans(img,ioff,joff,2*uv+i,j);
          
          for(ii=0;ii<4;ii++)
            for(jj=0;jj<4;jj++)
            {
              dec_picture->imgUV[uv][j4+jj][i4+ii]=img->m7[ii][jj];
            }
        }
      }
    }
  }

  return 0;
}


