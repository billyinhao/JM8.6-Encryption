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
 * \file image.c
 *
 * \brief
 *    Decode a Slice
 *
 * \author
 *    Main contributors (see contributors.h for copyright, address and affiliation details)
 *    - Inge Lille-Lang鴜               <inge.lille-langoy@telenor.com>
 *    - Rickard Sjoberg                 <rickard.sjoberg@era.ericsson.se>
 *    - Jani Lainema                    <jani.lainema@nokia.com>
 *    - Sebastian Purreiter             <sebastian.purreiter@mch.siemens.de>
 *    - Byeong-Moon Jeon                <jeonbm@lge.com>
 *    - Thomas Wedi                     <wedi@tnt.uni-hannover.de>
 *    - Gabi Blaettermann               <blaetter@hhi.de>
 *    - Ye-Kui Wang                     <wyk@ieee.org>
 *    - Antti Hallapuro                 <antti.hallapuro@nokia.com>
 *    - Alexis Tourapis                 <alexismt@ieee.org>
 ***********************************************************************
 */

#include "contributors.h"

#include <math.h>
#include <limits.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>


#include "global.h"
#include "errorconcealment.h"
#include "image.h"
#include "mbuffer.h"
#include "fmo.h"
#include "nalu.h"
#include "parsetcommon.h"
#include "parset.h"
#include "header.h"
#include "rtp.h"
#include "sei.h"
#include "output.h"
#include "biaridecod.h"
#include "mb_access.h"
#include "annexb.h"

#include "context_ini.h"
#include "cabac.h"
#include "loopfilter.h"

#include "vlc.h"

#include "erc_api.h"
extern objectBuffer_t *erc_object_list;
extern ercVariables_t *erc_errorVar;
extern frame erc_recfr;
extern int erc_mvperMB;
extern struct img_par *erc_img;

//extern FILE *p_out2;

extern StorablePicture **listX[6];
extern ColocatedParams *Co_located;

StorablePicture *dec_picture;

OldSliceParams old_slice;
void MbAffPostProc()
{
  byte temp[16][32];

  byte ** imgY  = dec_picture->imgY;
  byte ***imgUV = dec_picture->imgUV;

  int i, x, y, x0, y0, uv;
  for (i=0; i<(int)dec_picture->PicSizeInMbs; i+=2)
  {
    if (dec_picture->mb_field[i])
    {
      get_mb_pos(i, &x0, &y0);
      for (y=0; y<(2*MB_BLOCK_SIZE);y++)
        for (x=0; x<MB_BLOCK_SIZE; x++)
          temp[x][y] = imgY[y0+y][x0+x];

      for (y=0; y<MB_BLOCK_SIZE;y++)
        for (x=0; x<MB_BLOCK_SIZE; x++)
        {
          imgY[y0+(2*y)][x0+x]   = temp[x][y];
          imgY[y0+(2*y+1)][x0+x] = temp[x][y+MB_BLOCK_SIZE];
        }

      x0 = x0/2;
      y0 = y0/2;

      for (uv=0; uv<2; uv++)
      {
        for (y=0; y<(2*MB_BLOCK_SIZE/2);y++)
          for (x=0; x<MB_BLOCK_SIZE/2; x++)
            temp[x][y] = imgUV[uv][y0+y][x0+x];
          
        for (y=0; y<MB_BLOCK_SIZE/2;y++)
          for (x=0; x<MB_BLOCK_SIZE/2; x++)
          {
            imgUV[uv][y0+(2*y)][x0+x]   = temp[x][y];
            imgUV[uv][y0+(2*y+1)][x0+x] = temp[x][y+MB_BLOCK_SIZE/2];
          }
      }
    }
  }
}

/*!
 ***********************************************************************
 * \brief
 *    decodes one I- or P-frame
 *
 ***********************************************************************
 */

int decode_one_frame(struct img_par *img,struct inp_par *inp, struct snr_par *snr)
{
  int current_header;
  Slice *currSlice = img->currentSlice;

  img->current_slice_nr = 0;
  img->current_mb_nr = -4711;     // initialized to an impossible value for debugging -- correct value is taken from slice header
  currSlice->next_header = -8888; // initialized to an impossible value for debugging -- correct value is taken from slice header
  img->num_dec_mb = 0;
  img->newframe = 1;

  while ((currSlice->next_header != EOS && currSlice->next_header != SOP))
  {
    current_header = read_new_slice();//readesnew_slieeo对码流中的片数据进行读取(下面的片解码函数
   //decode_slice就是对读进的片数据进行解码)read-new-slice在此处主要作用是读取 NALU 单元数据 ，附
	//带做一些解码的准备工作
    if (current_header == EOS)
    {
      exit_picture();
      return EOS;
    }

    decode_slice(img, inp, current_header); //对之前读入的片数据进行解码

    img->newframe = 0;
    img->current_slice_nr++;
  }

  exit_picture();

  return (SOP);
}


/*!
 ************************************************************************
 * \brief
 *    Find PSNR for all three components.Compare decoded frame with
 *    the original sequence. Read inp->jumpd frames to reflect frame skipping.
 ************************************************************************
 */
void find_snr(
  struct snr_par  *snr,   //!< pointer to snr parameters
  StorablePicture *p,     //!< picture to be compared
  FILE *p_ref)            //!< open reference YUV file
{
  int i,j;
  int diff_y,diff_u,diff_v;
  int uv;
  int  status;

  // calculate frame number
  int  psnrPOC = active_sps->mb_adaptive_frame_field_flag ? p->poc /(input->poc_scale) : p->poc/(3-input->poc_scale);
//  int  psnrPOC = p->MbaffFrameFlag ? p->poc /(input->poc_scale) : p->poc/(3-input->poc_scale);

  // KS: Code below might work better if you have fields and a large (>1) poc offset between them
//  int  poc_diff=max(1,(p->bottom_poc - p->top_poc));
//  int  psnrPOC = active_sps->mb_adaptive_frame_field_flag ? p->poc /(input->poc_scale*poc_diff) : p->poc/((3-input->poc_scale)*poc_diff);

  if (psnrPOC==0 && img->psnr_number)
    img->idr_psnr_number=img->psnr_number + 1;
  img->psnr_number=max(img->psnr_number,img->idr_psnr_number+psnrPOC);
  
  frame_no = img->idr_psnr_number+psnrPOC;


  rewind(p_ref);

  for (i=0; i<frame_no; i++)
  {
    status = fseek (p_ref, (long) p->size_y* (long) (p->size_x*3/2), SEEK_CUR);
    if (status != 0)
    {
      snprintf(errortext, ET_SIZE, "Error in seeking frame number: %d", frame_no);
      fprintf(stderr, "%s", errortext);
      return;
//    snprintf(errortext, ET_SIZE, "Error in seeking frame number: %d", frame_no);
//    error(errortext, 500);
    }
  }

  for (j=0; j < p->size_y; j++)
    for (i=0; i < p->size_x; i++)
      imgY_ref[j][i]=fgetc(p_ref);

  for (uv=0; uv < 2; uv++)
    for (j=0; j < p->size_y_cr ; j++)
      for (i=0; i < p->size_x_cr; i++)
        imgUV_ref[uv][j][i]=fgetc(p_ref);

  img->quad[0]=0;
  diff_y=0;
  for (j=0; j < p->size_y; ++j)
  {
    for (i=0; i < p->size_x; ++i)
    {
      diff_y += img->quad[abs(p->imgY[j][i]-imgY_ref[j][i])];
    }
  }

  // Chroma
  diff_u=0;
  diff_v=0;

  for (j=0; j < p->size_y_cr; ++j)
  {
    for (i=0; i < p->size_x_cr; ++i)
    {
      diff_u += img->quad[abs(imgUV_ref[0][j][i]-p->imgUV[0][j][i])];
      diff_v += img->quad[abs(imgUV_ref[1][j][i]-p->imgUV[1][j][i])];
    }
  }

/*  if (diff_y == 0)
      diff_y = 1;
  if (diff_u == 0)
      diff_u = 1;
  if (diff_v == 0)
      diff_v = 1; */

  // Collecting SNR statistics
  if (diff_y != 0)
    snr->snr_y=(float)(10*log10(65025*(float)(p->size_x)*(p->size_y)/(float)diff_y));        // luma snr for current frame
  else
    snr->snr_y=0;
  if (diff_u != 0)
    snr->snr_u=(float)(10*log10(65025*(float)(p->size_x)*(p->size_y)/(float)(4*diff_u)));    //  chroma snr for current frame
  else
    snr->snr_u=0;
  if (diff_v != 0)
    snr->snr_v=(float)(10*log10(65025*(float)(p->size_x)*(p->size_y)/(float)(4*diff_v)));    //  chroma snr for current frame
  else
    snr->snr_v=0;

  if (img->number == 0) // first
  {
    snr->snr_ya=snr->snr_y1=snr->snr_y;                                                        // keep luma snr for first frame
    snr->snr_ua=snr->snr_u1=snr->snr_u;                                                        // keep chroma snr for first frame
    snr->snr_va=snr->snr_v1=snr->snr_v;                                                        // keep chroma snr for first frame
  
  }
  else
  {
    snr->snr_ya=(float)(snr->snr_ya*(img->number+Bframe_ctr)+snr->snr_y)/(img->number+Bframe_ctr+1); // average snr chroma for all frames
    snr->snr_ua=(float)(snr->snr_ua*(img->number+Bframe_ctr)+snr->snr_u)/(img->number+Bframe_ctr+1); // average snr luma for all frames
    snr->snr_va=(float)(snr->snr_va*(img->number+Bframe_ctr)+snr->snr_v)/(img->number+Bframe_ctr+1); // average snr luma for all frames
  } 
}


/*!
 ************************************************************************
 * \brief
 *    Interpolation of 1/4 subpixel
 ************************************************************************
 */
//++ 参见标准 8.4.2.2.1
//**如果是帧间编码，则还要判断是P或是B帧。找到参考帧的匹配块后，用get_block()像素内插恢复像素值。*****
void get_block(int ref_frame, StorablePicture **list, int x_pos, int y_pos, struct img_par *img, int block[BLOCK_SIZE][BLOCK_SIZE])
{

  int dx, dy;
  int x, y;
  int i, j;
  int maxold_x,maxold_y;
  int result;
  int pres_x;
  int pres_y; 
  int tmp_res[4][9];
  static const int COEF[6] = {    1, -5, 20, 20, -5, 1  };

  dx = x_pos&3;	//++ 当前4*4块的参考块的左上角像素（不一定在整像素位置）与其左方最邻近整像素的距离（以1/4像素距离为单位）
  dy = y_pos&3;	//++ 当前4*4块的参考块的左上角像素（不一定在整像素位置）与其上方最邻近整像素的距离（以1/4像素距离为单位）
  x_pos = (x_pos-dx)/4;	//++ 当前4*4块的参考块的左上角整像素的左方最邻近整像素在参考图像帧中的横坐标（以整像素距离为单位）
  y_pos = (y_pos-dy)/4;	//++ 当前4*4块的参考块的左上角整像素的上方最邻近整像素在参考图像帧中的纵坐标（以整像素距离为单位）

  maxold_x = dec_picture->size_x-1;
  maxold_y = dec_picture->size_y-1;

  if (dec_picture->mb_field[img->current_mb_nr])
    maxold_y = dec_picture->size_y/2 - 1;

  if (dx == 0 && dy == 0) {  /* fullpel position */	//++ 参考像素横纵坐标都在整像素点位置
    for (j = 0; j < BLOCK_SIZE; j++)
      for (i = 0; i < BLOCK_SIZE; i++)
        block[i][j] = list[ref_frame]->imgY[max(0,min(maxold_y,y_pos+j))][max(0,min(maxold_x,x_pos+i))];
  }
  else { /* other positions */

    if (dy == 0) { /* No vertical interpolation */	//++ 参考像素纵坐标在整像素点位置
	  //++ 水平方向6抽头滤波（1/2精度插值）
      for (j = 0; j < BLOCK_SIZE; j++) {
        for (i = 0; i < BLOCK_SIZE; i++) {
			for (result = 0, x = -2; x < 4; x++)
			  result += list[ref_frame]->imgY[max(0,min(maxold_y,y_pos+j))][max(0,min(maxold_x,x_pos+i+x))]*COEF[x+2];
		  block[i][j] = max(0, min(255, (result+16)/32));	//++ 四舍五入并限幅
		}
	  }

      if ((dx&1) == 1) {
		//++ 如果参考像素横坐标在1/4像素点位置，则进行水平方向1/4精度插值
        for (j = 0; j < BLOCK_SIZE; j++)																									//++++++++++++++++++++++++++
          for (i = 0; i < BLOCK_SIZE; i++)																									//++
            block[i][j] = (block[i][j] + list[ref_frame]->imgY[max(0,min(maxold_y,y_pos+j))][max(0,min(maxold_x,x_pos+i+dx/2))] +1 )/2;		//++
      }																																		//++++++++++++++++++++++++++
    }
    else if (dx == 0) {  /* No horizontal interpolation */	//++ 参考像素横坐标在整像素点位置
	  //++ 垂直方向6抽头滤波（1/2精度插值）
      for (j = 0; j < BLOCK_SIZE; j++) {
        for (i = 0; i < BLOCK_SIZE; i++) {
          for (result = 0, y = -2; y < 4; y++)
            result += list[ref_frame]->imgY[max(0,min(maxold_y,y_pos+j+y))][max(0,min(maxold_x,x_pos+i))]*COEF[y+2];
          block[i][j] = max(0, min(255, (result+16)/32));	//++ 四舍五入并限幅
        }
      }

      if ((dy&1) == 1) {
		//++ 如果参考像素纵坐标在1/4像素点位置，则进行垂直方向1/4精度插值
        for (j = 0; j < BLOCK_SIZE; j++)
          for (i = 0; i < BLOCK_SIZE; i++)
           block[i][j] = (block[i][j] + list[ref_frame]->imgY[max(0,min(maxold_y,y_pos+j+dy/2))][max(0,min(maxold_x,x_pos+i))] +1 )/2;
      }
    }
    else if (dx == 2) {  /* Vertical & horizontal interpolation */	//++ 参考像素横坐标在1/2像素点位置
	  //++ 水平方向6抽头滤波（1/2精度插值）
      for (j = -2; j < BLOCK_SIZE+3; j++) {
        for (i = 0; i < BLOCK_SIZE; i++)
          for (tmp_res[i][j+2] = 0, x = -2; x < 4; x++)
            tmp_res[i][j+2] += list[ref_frame]->imgY[max(0,min(maxold_y,y_pos+j))][max(0,min(maxold_x,x_pos+i+x))]*COEF[x+2];
      }
	  //++ 垂直方向6抽头滤波（1/2精度插值）
      for (j = 0; j < BLOCK_SIZE; j++) {
        for (i = 0; i < BLOCK_SIZE; i++) {
          for (result = 0, y = -2; y < 4; y++)
            result += tmp_res[i][j+y+2]*COEF[y+2];
          block[i][j] = max(0, min(255, (result+512)/1024));	//++ 四舍五入并限幅
        } 
      }

      if ((dy&1) == 1) {
		//++ 如果参考像素纵坐标在1/4像素点位置，则进行垂直方向1/4精度插值
        for (j = 0; j < BLOCK_SIZE; j++)
          for (i = 0; i < BLOCK_SIZE; i++)
            block[i][j] = (block[i][j] + max(0, min(255, (tmp_res[i][j+2+dy/2]+16)/32)) +1 )/2;
      }
    }
    else if (dy == 2) {  /* Horizontal & vertical interpolation */	//++ 参考像素纵坐标在1/2像素点位置
	  //++ 垂直方向6抽头滤波（1/2精度插值）
      for (j = 0; j < BLOCK_SIZE; j++) {
        for (i = -2; i < BLOCK_SIZE+3; i++)
          for (tmp_res[j][i+2] = 0, y = -2; y < 4; y++)
            tmp_res[j][i+2] += list[ref_frame]->imgY[max(0,min(maxold_y,y_pos+j+y))][max(0,min(maxold_x,x_pos+i))]*COEF[y+2];
      }
	  //++ 水平方向6抽头滤波（1/2精度插值）
      for (j = 0; j < BLOCK_SIZE; j++) {
        for (i = 0; i < BLOCK_SIZE; i++) {
          for (result = 0, x = -2; x < 4; x++)
            result += tmp_res[j][i+x+2]*COEF[x+2];
          block[i][j] = max(0, min(255, (result+512)/1024));	//++ 四舍五入并限幅
        }
      }

      if ((dx&1) == 1) {
		//++ 如果参考像素横坐标在1/4像素点位置，则进行水平方向1/4精度插值
        for (j = 0; j < BLOCK_SIZE; j++)
          for (i = 0; i < BLOCK_SIZE; i++)
            block[i][j] = (block[i][j] + max(0, min(255, (tmp_res[j][i+2+dx/2]+16)/32))+1)/2;
      }
    }
    else {  /* Diagonal interpolation */	//++ 参考像素横纵坐标都在1/4像素点位置

      for (j = 0; j < BLOCK_SIZE; j++) {
        for (i = 0; i < BLOCK_SIZE; i++) {
          pres_y = dy == 1 ? y_pos+j : y_pos+j+1;
          pres_y = max(0,min(maxold_y,pres_y));
          for (result = 0, x = -2; x < 4; x++)
            result += list[ref_frame]->imgY[pres_y][max(0,min(maxold_x,x_pos+i+x))]*COEF[x+2];
          block[i][j] = max(0, min(255, (result+16)/32));	//++ 四舍五入并限幅
        }
      }

      for (j = 0; j < BLOCK_SIZE; j++) {
        for (i = 0; i < BLOCK_SIZE; i++) {
          pres_x = dx == 1 ? x_pos+i : x_pos+i+1;
          pres_x = max(0,min(maxold_x,pres_x));
          for (result = 0, y = -2; y < 4; y++)
            result += list[ref_frame]->imgY[max(0,min(maxold_y,y_pos+j+y))][pres_x]*COEF[y+2];
          block[i][j] = (block[i][j] + max(0, min(255, (result+16)/32)) +1 ) / 2;	//++ 四舍五入并限幅
        }
      }

    }
  }

}


void reorder_lists(int currSliceType, Slice * currSlice)
{

  if ((currSliceType != I_SLICE)&&(currSliceType != SI_SLICE))
  {
    if (currSlice->ref_pic_list_reordering_flag_l0)
    {
      reorder_ref_pic_list(listX[0], &listXsize[0], 
                           img->num_ref_idx_l0_active - 1, 
                           currSlice->remapping_of_pic_nums_idc_l0, 
                           currSlice->abs_diff_pic_num_minus1_l0, 
                           currSlice->long_term_pic_idx_l0);
    }
    if (NULL == listX[0][img->num_ref_idx_l0_active-1])
    {
      error("number of entries in list 0 smaller than num_ref_idx_l0_active_minus1",500);
    }
    // that's a definition
    listXsize[0] = img->num_ref_idx_l0_active;
  }
  if (currSliceType == B_SLICE)
  {
    if (currSlice->ref_pic_list_reordering_flag_l1)
    {
      reorder_ref_pic_list(listX[1], &listXsize[1], 
                           img->num_ref_idx_l1_active - 1, 
                           currSlice->remapping_of_pic_nums_idc_l1, 
                           currSlice->abs_diff_pic_num_minus1_l1, 
                           currSlice->long_term_pic_idx_l1);
    }
    if (NULL == listX[1][img->num_ref_idx_l1_active-1])
    {
      error("number of entries in list 1 smaller than num_ref_idx_l1_active_minus1",500);
    }
    // that's a definition
    listXsize[1] = img->num_ref_idx_l1_active;
  }

  free_ref_pic_list_reordering_buffer(currSlice);
}


/*!
 ************************************************************************
 * \brief
 *    initialize ref_pic_num array
 ************************************************************************
 */
void set_ref_pic_num()  //set_ref_pic_num()的作用在于初始化参考帧序列
{
  int i,j;
  
  int slice_id=img->current_slice_nr;

  for (i=0;i<listXsize[LIST_0];i++)
  {
    dec_picture->ref_pic_num        [slice_id][LIST_0][i]=listX[LIST_0][i]->poc * 2 + ((listX[LIST_0][i]->structure==BOTTOM_FIELD)?1:0) ;	//??? 为什么*2
    dec_picture->frm_ref_pic_num    [slice_id][LIST_0][i]=listX[LIST_0][i]->frame_poc * 2; 
    dec_picture->top_ref_pic_num    [slice_id][LIST_0][i]=listX[LIST_0][i]->top_poc * 2; 
    dec_picture->bottom_ref_pic_num [slice_id][LIST_0][i]=listX[LIST_0][i]->bottom_poc * 2 + 1; 
    //printf("POCS %d %d %d %d ",listX[LIST_0][i]->frame_poc,listX[LIST_0][i]->bottom_poc,listX[LIST_0][i]->top_poc,listX[LIST_0][i]->poc);
    //printf("refid %d %d %d %d\n",(int) dec_picture->frm_ref_pic_num[LIST_0][i],(int) dec_picture->top_ref_pic_num[LIST_0][i],(int) dec_picture->bottom_ref_pic_num[LIST_0][i],(int) dec_picture->ref_pic_num[LIST_0][i]);
  }

  for (i=0;i<listXsize[LIST_1];i++)
  {
    dec_picture->ref_pic_num        [slice_id][LIST_1][i]=listX[LIST_1][i]->poc  *2 + ((listX[LIST_1][i]->structure==BOTTOM_FIELD)?1:0);
    dec_picture->frm_ref_pic_num    [slice_id][LIST_1][i]=listX[LIST_1][i]->frame_poc * 2; 
    dec_picture->top_ref_pic_num    [slice_id][LIST_1][i]=listX[LIST_1][i]->top_poc * 2; 
    dec_picture->bottom_ref_pic_num [slice_id][LIST_1][i]=listX[LIST_1][i]->bottom_poc * 2 + 1; 
  }

  if (img->structure==FRAME)
    for (j=2;j<6;j++)
      for (i=0;i<listXsize[j];i++)
      {    
        dec_picture->ref_pic_num        [slice_id][j][i] = listX[j][i]->poc * 2 + ((listX[j][i]->structure==BOTTOM_FIELD)?1:0);
        dec_picture->frm_ref_pic_num    [slice_id][j][i] = listX[j][i]->frame_poc * 2 ;
        dec_picture->top_ref_pic_num    [slice_id][j][i] = listX[j][i]->top_poc * 2 ;
        dec_picture->bottom_ref_pic_num [slice_id][j][i] = listX[j][i]->bottom_poc * 2 + 1;
      }

}


/*!
 ************************************************************************
 * \brief
 *    Reads new slice from bit_stream
 ************************************************************************
 */
/****************************************************************************************
read_new_slice()过程读取编码片中的数据（注：JM模型解码的过程：在进行相应的初始工作之后
循环调用Decode--one--frame()函数，对视频编码序列的每一帧进行解码处理。2帧解码过程中，先通
过readesnew_slieeo对码流中的片数据进行读取，然后进行片解码decode一llceo。接着判断该片是否
解码完成，如果没有完成，则继续读循环取下一片数据，再解码。如果所有片解码完成，则对该图像进
行去方块效应滤波。滤波是对图像中每一个宏块进行的，即DeblockMbO函数。该函数又包含获取边界强
度GctstrengthO和滤波过程Edge功oPO。进行滤波之后，把已解码帧放入解码帧缓冲器(DPB)，作为后续
解码过程的参考帧。放入解码帧缓冲器的中的参考帧是按队列存放的。默认情况下，解码帧缓冲器最多
存放16帧图像。当解码帧缓充器中的参考帧存放满而没有空间时，如有新的帧要进入则输出最早进入队
列的帧，其它帧按顺序在队列中前移一位。
******************************************************************************************/

int read_new_slice()  //read_new_slice 函数的主要作用是读取 NALU 单元数据,附带做一些解码的准备工作  
{
  NALU_t *nalu = AllocNALU(MAX_CODED_FRAME_SIZE);//通过AllocNALU()为NAL单元分配一个大小为
                                                 //MAXeeCODEDesFRAME--SIzE的堆空间，用来存放己编码的帧

  int current_header;
  int ret;
  int BitsUsedByHeader;
  Slice *currSlice = img->currentSlice;
  Bitstream *currStream;

  int slice_id_a, slice_id_b, slice_id_c;
  int redundant_pic_cnt_b, redundant_pic_cnt_c;
  long ftell_position, expected_slice_type;
  
//  int i;
  expected_slice_type = NALU_TYPE_DPA;

  while (1)
  {
    ftell_position = ftell(bits);	//++ ftell()得到某文件指针的当前位置

    if (input->FileFormat == PAR_OF_ANNEXB)
      ret=GetAnnexbNALU (nalu);	//++ GetAnnexbNALU (nalu)返回的值是正在使用的NAL单元的长度（包括开始前缀码），字节码流
    else
      ret=GetRTPNALU (nalu);//GetRTPNALU 处理 RTP 格式码流
	//字节流格式的码流主要用于存储，例如制作 DVD（当然现在的 DVD 还不是用 H.264),RTP 
	//格式码流主要用于网络传送，例如在线看电影 

    //我们知道H.264 码流的第一个 NALU 是 SPS ,所以 7.3 部分接下来就描述了 SPS 的语法结构,
	//解码器下一步工作当然是要解码 SPS(SPS 的结构，请看 7.3.2.1 小节,加粗的字体就是被编码的语法元素，
	//也就是解码器必须解码的语法元素 ),解码器必须按照 7.3 部分中加粗字体的顺序来解码 ,
	//对于 SPS，按照 7.3.2.1 的规定，解码器必须首先解码 profile_idc 

    //In some cases, zero_byte shall be present. If current NALU is a VCL NALU, we can't tell
    //whether it is the first VCL NALU at this point, so only non-VCL NAL unit is checked here.
    CheckZeroByteNonVCL(nalu, &ret);	//++ 该函数同时执行了对已解码的非VCL NAL单元计数

    NALUtoRBSP(nalu);	//++ 剔除比特流中的仿校验字节（0X03）
//    printf ("nalu->len %d\n", nalu->len);
    
    if (ret < 0)
      printf ("Error while getting the NALU in file format %s, exit\n", input->FileFormat==PAR_OF_ANNEXB?"Annex B":"RTP");
    if (ret == 0)
    {
//      printf ("read_new_slice: returning %s\n", "EOS");
      if(expected_slice_type != NALU_TYPE_DPA)
      {
        /* oops... we found the next slice, go back! */
        fseek(bits, ftell_position, SEEK_SET);
        FreeNALU(nalu);
        return current_header;
      }
      else
        return EOS;
    }

    // Got a NALU
    if (nalu->forbidden_bit)
    {
      printf ("Found NALU w/ forbidden_bit set, bit error?  Let's try...\n");
    }

    switch (nalu->nal_unit_type)//判断NAL单元类型，根据nal一nit--type的不同取值，进入不同的解码方式。
		                        //接着是一些初始化，然后deeod_poc()函数计算图像序列号 Poe(Pietureorder
                                //Count)，根据pic--order--cnt--type值的不同，会采用不同的算法来计算Poe

    {
      case NALU_TYPE_SLICE: //case语句，根据 NALU 的头解析出来的 nal_unit_type,
		                   //解码器会进行不同 NALU 单元的处理,这个 switch 语句中包括对 一般 slice 的处理 ,
		                    //对 IDR 的处理,对 DP1、DPB、DPC 的处理,包括对 SPS、PPS 的处理等等 
		  //解码器现在遇到的第一个 NALU 是 SPS,因此解码器经过解析之后的 nal_unit_type 的值必然为 7,
		  //解码器就会跳转到 case NALU_TYPE_SPS并进入 ProcessSPS 函数 

		  //我们知道 H.264 码流的第一个 NALU 是 SPS,而H.264 码流的第二个 NALU 是是 PPS，所以在处理完了
		  //case NALU_TYPE_SPS相关的操作后，解码器下一步必然会跳转到 read_new_slice 函数中的 case NALU_TYPE_PPS,
		  //并调用 ProcessPPS 函数,之后进这个函数看一下

      case NALU_TYPE_IDR:
		  //H.264 码流的第三个 NALU 假设是 IDR（因为有可能是定界符，看一下毕厚杰的书），所以解码器在处理完前两个NALU
		  //后，就要跳转到 read_new_slice 函数中的 case NALU_TYPE_IDR 
        img->idr_flag = (nalu->nal_unit_type == NALU_TYPE_IDR);
        img->nal_reference_idc = nalu->nal_reference_idc;
        img->disposable_flag = (nalu->nal_reference_idc == NALU_PRIORITY_DISPOSABLE);
        currSlice->dp_mode = PAR_DP_1;	//++ dp_mode：数据分割模式；PAR_DP_1=0：没有数据分割
        currSlice->max_part_nr = 1;
        currSlice->ei_flag = 0;	//++ 该处赋值直接影响decode_slice()函数中对decode_one_slice()函数的调用
									//++ 该值不为0，表明当前片出错，解码程序将忽略当前片的解码过程，而使用错误隐藏
        currStream = currSlice->partArr[0].bitstream;
        currStream->ei_flag = 0;	//++ 此处的赋值为最终赋值，以后不再改变。该值将对每个宏块的ei_flag产生影响
									//++ 参见macroblock.c文件read_one_macroblock()函数的如下语句：
									//++		:if(!dP->bitstream->ei_flag)		:currMB->ei_flag = 0;
									//++ 该值还在macroblock.c文件if(IS_INTRA (currMB) && dP->bitstream->ei_flag && img->number)中用到
        currStream->frame_bitoffset = currStream->read_len = 0;
        memcpy (currStream->streamBuffer, &nalu->buf[1], nalu->len-1);
        currStream->code_len = currStream->bitstream_length = RBSPtoSODB(currStream->streamBuffer, nalu->len-1);

        // Some syntax of the Slice Header depends on the parameter set, which depends on
        // the parameter set ID of the SLice header.  Hence, read the pic_parameter_set_id
        // of the slice header first, then setup the active parameter sets, and then read
        // the rest of the slice header

/*
因为 H.264 中的图像都是按照片来组织的,因此片作为最大的语法结构而存在,不考虑片组,
H.264 码流是一个个 NALU 拼在一起的,首尾相连 , 每个 NALU 由头和体组成,NALU 头就是
刚才我们看的 7.3.1 小节中的三个粗体字forbidden_zero_bit（1比特）、nal_ref_idc（2比特）、
nal_unit_type（5比特）, NALU 体就是 RBSP 数据[EBSP：在RBSP的基础上增加了防止伪起始码字
节(仿校验字节)（0X03)],而EBSP 数据又是由一个个片组成的，片又由片头和片数据组成，
所以解码器取出整个片数据后，接下来必然要处理片头. 我们已经进入到了 FirstPartOfSliceHeader 
函数，从这个函数的名字上，我们就已经知道它在做什么了.那么解码器在这个函数里具体该怎么做呢？ 
 标准 7.3.3 小节已经做了规定.所以解码器必然是按照标准 7.3.3 小节规定的顺序来做的,标准 7.3.3 
小节规定的其他内容也都是在 RestOfSliceHeader 函数里完成，虽然标准 7.3.3 小节规定的内容被分到
了两个函数中进行处理,但是解码器必然不可能打乱标准 7.3.3 小节规定的顺序,所以虽然在解码器代码
中的 FirstPartOfSliceHeader 函数与 RestOfSliceHeader 函数之间插入了一个 UseParameterSet 函数,  
但是 UseParameterSet 函数绝对不可能读码流,因为 FirstPartOfSliceHeader 函数处理的最后一个语法
元素是 pic_parameter_set_id ，RestOfSliceHeader 函数处理的第一个语法元素是 frame_num ,
而标准 7.3.3 小节规定语法元素 pic_parameter_set_id 之后紧接着是语法元素frame_num ，所以
如果 UseParameterSet 函数读了码流就必然导致比特错误
*/
        BitsUsedByHeader = FirstPartOfSliceHeader();	//++ 参见标准7.3.3，进入FirstPartOfSliceHeader()
        UseParameterSet (currSlice->pic_parameter_set_id);//UseParameterSet 函数做了什么呢？
		//进去之后可以看见它主要是激活 SPS 和 PPS,因为后面语法元素的解码会使用到 SPS、PPS 里的参数，
		//所以在解码后续语法元素之前必须先找到正确的 SPS、PPS,请看标准 7.3.3 中的 ref_pic_list_reordering( ) 
		//7.3.3 规定了解码完片头的语法元素，要进行 ref_pic_list_reordering 操作，那么解码器也应该这样做
        BitsUsedByHeader+= RestOfSliceHeader ();	//++ 参见标准7.3.3，进入RestOfSliceHeader ()
		//++ BitsUsedByHeader在程序中没有实际用处，而且BitsUsedByHeader+= RestOfSliceHeader ()
		//++ 重复计算了FirstPartOfSliceHeader()所用到的比特数。因为在FirstPartOfSliceHeader()
		//++ 之后，变量UsedBits值并未被置零就代入RestOfSliceHeader()运算，从而RestOfSliceHeader ()
		//++ 在返回时，BitsUsedByHeader+= RestOfSliceHeader()多加了一个BitsUsedByHeader值

        FmoInit (active_pps, active_sps);

        if(is_new_picture())
        {
          init_picture(img, input);
          
          current_header = SOP;
          //check zero_byte if it is also the first NAL unit in the access unit
          CheckZeroByteVCL(nalu, &ret);
        }
        else
          current_header = SOS;
  
        init_lists(img->type, img->currentSlice->structure);
        reorder_lists (img->type, img->currentSlice);

        if (img->structure==FRAME)
        {
          init_mbaff_lists();
        }

/*        if (img->frame_num==1) // write a reference list
        {
          count ++;
          if (count==1)
            for (i=0; i<listXsize[0]; i++)
              write_picture(listX[0][i], p_out2);
        }
*/

        // From here on, active_sps, active_pps and the slice header are valid
        if (img->MbaffFrameFlag)
          img->current_mb_nr = currSlice->start_mb_nr << 1;
        else
          img->current_mb_nr = currSlice->start_mb_nr;

        if (active_pps->entropy_coding_mode_flag)
        {
          int ByteStartPosition = currStream->frame_bitoffset/8;
          if (currStream->frame_bitoffset%8 != 0) 
          {
            ByteStartPosition++;
          }
          arideco_start_decoding (&currSlice->partArr[0].de_cabac, currStream->streamBuffer, ByteStartPosition, &currStream->read_len, img->type);
        }
// printf ("read_new_slice: returning %s\n", current_header == SOP?"SOP":"SOS");
        FreeNALU(nalu); //释放分配的NAL单元堆空间
        return current_header;
        break;
      case NALU_TYPE_DPA:
        //! The state machine here should follow the same ideas as the old readSliceRTP()
        //! basically:
        //! work on DPA (as above)
        //! read and process all following SEI/SPS/PPS/PD/Filler NALUs
        //! if next video NALU is dpB, 
        //!   then read and check whether it belongs to DPA, if yes, use it
        //! else
        //!   ;   // nothing
        //! read and process all following SEI/SPS/PPS/PD/Filler NALUs
        //! if next video NALU is dpC
        //!   then read and check whether it belongs to DPA (and DPB, if present), if yes, use it, done
        //! else
        //!   use the DPA (and the DPB if present)

        /* 
            LC: inserting the code related to DP processing, mainly copying some of the parts
            related to NALU_TYPE_SLICE, NALU_TYPE_IDR.
        */

        if(expected_slice_type != NALU_TYPE_DPA)
        {
          /* oops... we found the next slice, go back! */
          fseek(bits, ftell_position, SEEK_SET);
          FreeNALU(nalu);
          return current_header;
        }

        img->idr_flag          = (nalu->nal_unit_type == NALU_TYPE_IDR);
        img->nal_reference_idc = nalu->nal_reference_idc;
        img->disposable_flag   = (nalu->nal_reference_idc == NALU_PRIORITY_DISPOSABLE);
        currSlice->dp_mode     = PAR_DP_3;
        currSlice->max_part_nr = 3;
        currSlice->ei_flag     = 0;
        currStream             = currSlice->partArr[0].bitstream;
        currStream->ei_flag    = 0;
        currStream->frame_bitoffset = currStream->read_len = 0;
        memcpy (currStream->streamBuffer, &nalu->buf[1], nalu->len-1);
        currStream->code_len = currStream->bitstream_length = RBSPtoSODB(currStream->streamBuffer, nalu->len-1);	//++ 剔除停止比特和填充比特
        
        BitsUsedByHeader     = FirstPartOfSliceHeader();
        UseParameterSet (currSlice->pic_parameter_set_id);
        BitsUsedByHeader    += RestOfSliceHeader ();
        
        FmoInit (active_pps, active_sps);

        if(is_new_picture())
        {
          init_picture(img, input);
          current_header = SOP;
          CheckZeroByteVCL(nalu, &ret);
        }
        else
          current_header = SOS;

        
        init_lists(img->type, img->currentSlice->structure);
        reorder_lists (img->type, img->currentSlice);
        
        if (img->structure==FRAME)
        {
          init_mbaff_lists();
        }

        // From here on, active_sps, active_pps and the slice header are valid
        if (img->MbaffFrameFlag)
          img->current_mb_nr = currSlice->start_mb_nr << 1;
        else
          img->current_mb_nr = currSlice->start_mb_nr;


        /* 
           LC:
              Now I need to read the slice ID, which depends on the value of 
              redundant_pic_cnt_present_flag (pag.49). 
        */

        slice_id_a  = ue_v("NALU:SLICE_A slice_idr", currStream);
        if (active_pps->entropy_coding_mode_flag)
        {
          int ByteStartPosition = currStream->frame_bitoffset/8;
          if (currStream->frame_bitoffset%8 != 0) 
          {
            ByteStartPosition++;
          }
          arideco_start_decoding (&currSlice->partArr[0].de_cabac, currStream->streamBuffer, ByteStartPosition, &currStream->read_len, img->type);
        }
// printf ("read_new_slice: returning %s\n", current_header == SOP?"SOP":"SOS");
        break;
      case NALU_TYPE_DPB:
        /* LC: inserting the code related to DP processing */

        currStream             = currSlice->partArr[1].bitstream;
        currStream->ei_flag    = 0;
        currStream->frame_bitoffset = currStream->read_len = 0;
        memcpy (currStream->streamBuffer, &nalu->buf[1], nalu->len-1);
        currStream->code_len = currStream->bitstream_length = RBSPtoSODB(currStream->streamBuffer, nalu->len-1);

        slice_id_b  = ue_v("NALU:SLICE_B slice_idr", currStream);
        if (active_pps->redundant_pic_cnt_present_flag)
          redundant_pic_cnt_b = ue_v("NALU:SLICE_B redudand_pic_cnt", currStream);
        else
          redundant_pic_cnt_b = 0;
        
        /*  LC: Initializing CABAC for the current data stream. */

        if (active_pps->entropy_coding_mode_flag)
        {
          int ByteStartPosition = currStream->frame_bitoffset/8;
          if (currStream->frame_bitoffset % 8 != 0) 
            ByteStartPosition++;
          
          arideco_start_decoding (&currSlice->partArr[1].de_cabac, currStream->streamBuffer, 
            ByteStartPosition, &currStream->read_len, img->type);
          
        }

        /* LC: resilience code to be inserted */
        /*         FreeNALU(nalu); */
        /*         return current_header; */

        break;
      case NALU_TYPE_DPC:
        /* LC: inserting the code related to DP processing */
        currStream             = currSlice->partArr[2].bitstream;
        currStream->ei_flag    = 0;
        currStream->frame_bitoffset = currStream->read_len = 0;
        memcpy (currStream->streamBuffer, &nalu->buf[1], nalu->len-1);
        currStream->code_len = currStream->bitstream_length = RBSPtoSODB(currStream->streamBuffer, nalu->len-1);
        
        slice_id_c  = ue_v("NALU:SLICE_C slice_idr", currStream);
        if (active_pps->redundant_pic_cnt_present_flag)
          redundant_pic_cnt_c = ue_v("NALU:SLICE_C redudand_pic_cnt", currStream);
        else
          redundant_pic_cnt_c = 0;
        
        /* LC: Initializing CABAC for the current data stream. */

        if (active_pps->entropy_coding_mode_flag)
        {
          int ByteStartPosition = currStream->frame_bitoffset/8;
          if (currStream->frame_bitoffset % 8 != 0) 
            ByteStartPosition++;
          
          arideco_start_decoding (&currSlice->partArr[2].de_cabac, currStream->streamBuffer, 
            ByteStartPosition, &currStream->read_len, img->type);
        }

        /* LC: resilience code to be inserted */

        FreeNALU(nalu);//释放分配的NAL单元堆空间，read--new_sliee过程结束。

        return current_header;

        break;
      case NALU_TYPE_SEI:
        printf ("read_new_slice: Found NALU_TYPE_SEI, len %d\n", nalu->len);
        InterpretSEIMessage(nalu->buf,nalu->len,img);
        break;
      case NALU_TYPE_PPS:
        ProcessPPS(nalu);
        break;

      case NALU_TYPE_SPS:
        ProcessSPS(nalu);
        break;
      case NALU_TYPE_AUD:
//        printf ("read_new_slice: Found 'Access Unit Delimiter' NAL unit, len %d, ignored\n", nalu->len);
        break;
      case NALU_TYPE_EOSEQ:
//        printf ("read_new_slice: Found 'End of Sequence' NAL unit, len %d, ignored\n", nalu->len);
        break;
      case NALU_TYPE_EOSTREAM:
//        printf ("read_new_slice: Found 'End of Stream' NAL unit, len %d, ignored\n", nalu->len);
        break;
      case NALU_TYPE_FILL:
        printf ("read_new_slice: Found NALU_TYPE_FILL, len %d\n", nalu->len);
        printf ("Skipping these filling bits, proceeding w/ next NALU\n");
        break;
      default:
        printf ("Found NALU type %d, len %d undefined, ignore NALU, moving on\n", nalu->nal_unit_type, nalu->len);
    }
  }
  FreeNALU(nalu);//最后,执行函数FreeNALUO，释放分配的NAL单元堆空间，read--new_sliee过程结束。

  return  current_header;
}


/*!
 ************************************************************************
 * \brief
 *    Initializes the parameters for a new picture
 ************************************************************************
 */
void init_picture(struct img_par *img, struct inp_par *inp)
{
  int i,k,l;
  if (dec_picture)
  {
    // this may only happen on slice loss
    exit_picture();
  }

  if (img->frame_num != img->pre_frame_num && img->frame_num != (img->pre_frame_num + 1) % img->MaxFrameNum) 
  {
    if (active_sps->gaps_in_frame_num_value_allowed_flag == 0)
    {
      /* Advanced Error Concealment would be called here to combat unintentional loss of pictures. */
      error("An unintentional loss of pictures occurs! Exit\n", 100);
    }
    fill_frame_num_gap(img);
  }
  img->pre_frame_num = img->frame_num;
  img->num_dec_mb = 0;

  //calculate POC
  decode_poc(img);//deeod_poc()函数计算图像序列号POC(picture order count),根据
                  //pic--order--cnt--type值的不同，会采用不同的算法来计算Poc
  //  dumppoc (img);

  if (img->structure==FRAME ||img->structure==TOP_FIELD)
  {
#ifdef WIN32
    _ftime (&(img->tstruct_start));             // start time ms
#else
    ftime (&(img->tstruct_start));              // start time ms
#endif
    time( &(img->ltime_start));                // start time s
  }

  dec_picture = alloc_storable_picture (img->structure, img->width, img->height, img->width_cr, img->height_cr);
  dec_picture->top_poc=img->toppoc;
  dec_picture->bottom_poc=img->bottompoc;
  dec_picture->frame_poc=img->framepoc;

  // reset all variables of the error concealment instance before decoding of every frame.
  // here the third parameter should, if perfectly, be equal to the number of slices per frame.
  // using little value is ok, the code will allocate more memory if the slice number is larger
  ercReset(erc_errorVar, img->PicSizeInMbs, img->PicSizeInMbs, dec_picture->size_x);
  erc_mvperMB = 0;

  switch (img->structure )
  {
  case TOP_FIELD:
    {
      dec_picture->poc=img->toppoc;
      img->number *= 2;
      break;
    }
  case BOTTOM_FIELD:
    {
      dec_picture->poc=img->bottompoc;
      img->number++;
      break;
    }
  case FRAME:
    {
      dec_picture->poc=img->framepoc;
      break;
    }
  default:
    error("img->structure not initialized", 235);
  }

  img->current_slice_nr=0;

  if (img->type > SI_SLICE)
  {
    set_ec_flag(SE_PTYPE);
    img->type = P_SLICE;  // concealed element
  }

  // CAVLC init
  for (i=0;i < (int)img->PicSizeInMbs; i++)
    for (k=0;k<4;k++)	//++ 代表每个8*8块分为四个4*4的块
      for (l=0;l<6;l++)	//++ 四个8*8亮度块+二个8*8色度块
        img->nz_coeff[i][k][l]=-1;  // CAVLC

  if(active_pps->constrained_intra_pred_flag)
  {
    for (i=0; i<(int)img->PicSizeInMbs; i++)
    {
      img->intra_block[i] = 1;
    }
  }

  // Set the slice_nr member of each MB to -1, to ensure correct when packet loss occurs
  // TO set Macroblock Map (mark all MBs as 'have to be concealed')
  for(i=0; i<(int)img->PicSizeInMbs; i++)
  {
    img->mb_data[i].slice_nr = -1; 
    img->mb_data[i].ei_flag = 1;
  }

  img->mb_y = img->mb_x = 0;
  img->block_y = img->pix_y = img->pix_c_y = 0; // define vertical positions
  img->block_x = img->pix_x = img->pix_c_x = 0; // define horizontal positions

  dec_picture->slice_type = img->type;
  dec_picture->used_for_reference = (img->nal_reference_idc != 0);
  dec_picture->idr_flag = img->idr_flag;
  dec_picture->no_output_of_prior_pics_flag = img->no_output_of_prior_pics_flag;
  dec_picture->long_term_reference_flag = img->long_term_reference_flag;
  dec_picture->adaptive_ref_pic_buffering_flag = img->adaptive_ref_pic_buffering_flag;

  dec_picture->dec_ref_pic_marking_buffer = img->dec_ref_pic_marking_buffer;
  img->dec_ref_pic_marking_buffer = NULL;

  dec_picture->MbaffFrameFlag = img->MbaffFrameFlag;
  dec_picture->PicWidthInMbs = img->PicWidthInMbs;
  dec_picture->pic_num = img->frame_num;
  dec_picture->frame_num = img->frame_num;
  dec_picture->coded_frame = (img->structure==FRAME);

  dec_picture->frame_mbs_only_flag = active_sps->frame_mbs_only_flag;
  dec_picture->frame_cropping_flag = active_sps->frame_cropping_flag;

  if (dec_picture->frame_cropping_flag)
  {
    dec_picture->frame_cropping_rect_left_offset   = active_sps->frame_cropping_rect_left_offset;
    dec_picture->frame_cropping_rect_right_offset  = active_sps->frame_cropping_rect_right_offset;
    dec_picture->frame_cropping_rect_top_offset    = active_sps->frame_cropping_rect_top_offset;
    dec_picture->frame_cropping_rect_bottom_offset = active_sps->frame_cropping_rect_bottom_offset;
  }

}

/*!
 ************************************************************************
 * \brief
 *    finish decoding of a picture, conceal errors and store it 
 *    into the DPB
 ************************************************************************
 */
void exit_picture()
{
  int ercStartMB;
  int ercSegment;
  frame recfr;
  unsigned int i;
  int structure, frame_poc, slice_type, refpic;

  int tmp_time;                   // time used by decoding the last frame

  // return if the last picture has already been finished
  if (dec_picture==NULL)
  {
    return;
  }

  //deblocking for frame or field

  DeblockPicture( img, dec_picture );

  if (dec_picture->MbaffFrameFlag)
    MbAffPostProc();

  recfr.yptr = &dec_picture->imgY[0][0];
  recfr.uptr = &dec_picture->imgUV[0][0][0];
  recfr.vptr = &dec_picture->imgUV[1][0][0];

  //! this is always true at the beginning of a picture
  ercStartMB = 0;
  ercSegment = 0;

  //! mark the start of the first segment
  if (!dec_picture->MbaffFrameFlag)
  {
    ercStartSegment(0, ercSegment, 0 , erc_errorVar);
    //! generate the segments according to the macroblock map
    for(i = 1; i<dec_picture->PicSizeInMbs; i++)
    {
      if(img->mb_data[i].ei_flag != img->mb_data[i-1].ei_flag)
      {
        ercStopSegment(i-1, ercSegment, 0, erc_errorVar); //! stop current segment
        
        //! mark current segment as lost or OK
        if(img->mb_data[i-1].ei_flag)
          ercMarkCurrSegmentLost(dec_picture->size_x, erc_errorVar);
        else
          ercMarkCurrSegmentOK(dec_picture->size_x, erc_errorVar);
        
        ercSegment++;  //! next segment
        ercStartSegment(i, ercSegment, 0 , erc_errorVar); //! start new segment
        ercStartMB = i;//! save start MB for this segment 
      }
    }
    //! mark end of the last segment
    ercStopSegment(dec_picture->PicSizeInMbs-1, ercSegment, 0, erc_errorVar);
    if(img->mb_data[i-1].ei_flag)
      ercMarkCurrSegmentLost(dec_picture->size_x, erc_errorVar);
    else
      ercMarkCurrSegmentOK(dec_picture->size_x, erc_errorVar);
    
    //! call the right error concealment function depending on the frame type.
    erc_mvperMB /= dec_picture->PicSizeInMbs;
    
    erc_img = img;
    if(dec_picture->slice_type == I_SLICE || dec_picture->slice_type == SI_SLICE) // I-frame
      ercConcealIntraFrame(&recfr, dec_picture->size_x, dec_picture->size_y, erc_errorVar);
    else
      ercConcealInterFrame(&recfr, erc_object_list, dec_picture->size_x, dec_picture->size_y, erc_errorVar);
  }

  if (img->structure == FRAME)         // buffer mgt. for frame mode
    frame_postprocessing(img, input);
  else
    field_postprocessing(img, input);   // reset all interlaced variables

  structure  = dec_picture->structure;
  slice_type = dec_picture->slice_type;
  frame_poc  = dec_picture->frame_poc;
  refpic     = dec_picture->used_for_reference;

  store_picture_in_dpb(dec_picture);
  dec_picture=NULL;

  if (img->last_has_mmco_5)
  {
    img->pre_frame_num = 0;
  }

  if ((structure==FRAME)||structure==BOTTOM_FIELD)
  {
    
#ifdef WIN32
    _ftime (&(img->tstruct_end));             // start time ms
#else
    ftime (&(img->tstruct_end));              // start time ms
#endif
    
    time( &(img->ltime_end));                // start time s

    tmp_time=(img->ltime_end*1000+img->tstruct_end.millitm) - (img->ltime_start*1000+img->tstruct_start.millitm);
    tot_time=tot_time + tmp_time;
    
    if(slice_type == I_SLICE) // I picture
      fprintf(stdout,"%3d(I)  %3d %5d %7.4f %7.4f %7.4f %5d\n",
      frame_no, frame_poc, img->qp,snr->snr_y,snr->snr_u,snr->snr_v,tmp_time);
    else if(slice_type == P_SLICE) // P pictures
      fprintf(stdout,"%3d(P)  %3d %5d %7.4f %7.4f %7.4f %5d\n",
      frame_no, frame_poc, img->qp,snr->snr_y,snr->snr_u,snr->snr_v,tmp_time);
    else if(slice_type == SP_SLICE) // SP pictures
      fprintf(stdout,"%3d(SP) %3d %5d %7.4f %7.4f %7.4f %5d\n",
      frame_no, frame_poc, img->qp,snr->snr_y,snr->snr_u,snr->snr_v,tmp_time);
    else if (slice_type == SI_SLICE)
      fprintf(stdout,"%3d(SI) %3d %5d %7.4f %7.4f %7.4f %5d\n",
      frame_no, frame_poc, img->qp,snr->snr_y,snr->snr_u,snr->snr_v,tmp_time);
    else if(refpic) // stored B pictures
      fprintf(stdout,"%3d(BS) %3d %5d %7.4f %7.4f %7.4f %5d\n",
      frame_no, frame_poc, img->qp,snr->snr_y,snr->snr_u,snr->snr_v,tmp_time);
    else // B pictures
      fprintf(stdout,"%3d(B)  %3d %5d %7.4f %7.4f %7.4f %5d\n",
      frame_no, frame_poc, img->qp,snr->snr_y,snr->snr_u,snr->snr_v,tmp_time);
    
    fflush(stdout);
    
    if(slice_type == I_SLICE || slice_type == SI_SLICE || slice_type == P_SLICE || refpic)   // I or P pictures
      img->number++;
    else
      Bframe_ctr++;    // B pictures
 
    g_nFrame++;
  }

  img->current_mb_nr = -4712;   // impossible value for debugging, StW
  img->current_slice_nr = 0;

}

/*!
 ************************************************************************
 * \brief
 *    write the encoding mode and motion vectors of current 
 *    MB to the buffer of the error concealment module.
 ************************************************************************
 */

void ercWriteMBMODEandMV(struct img_par *img,struct inp_par *inp)
{
  extern objectBuffer_t *erc_object_list;
  int i, ii, jj, currMBNum = img->current_mb_nr;
  int mbx = xPosMB(currMBNum,dec_picture->size_x), mby = yPosMB(currMBNum,dec_picture->size_x);	//++ 计算当前宏块的坐标
  objectBuffer_t *currRegion, *pRegion;
  Macroblock *currMB = &img->mb_data[currMBNum];
  int***  mv;

  currRegion = erc_object_list + (currMBNum<<2);

  if(img->type != B_SLICE) //non-B frame
  {
    for (i=0; i<4; i++)
    {
      pRegion             = currRegion + i;																//////////////////////////////////
      pRegion->regionMode = (currMB->mb_type  ==I16MB  ? REGMODE_INTRA      :								//++
                             currMB->b8mode[i]==IBLOCK ? REGMODE_INTRA_8x8  :								//++ 设置每个8*8块的预测模式
                             currMB->b8mode[i]==0      ? REGMODE_INTER_COPY :								//++
                             currMB->b8mode[i]==1      ? REGMODE_INTER_PRED : REGMODE_INTER_PRED_8x8);	//////////////////////////////////
      if (currMB->b8mode[i]==0 || currMB->b8mode[i]==IBLOCK)  // INTRA OR COPY
      {
        pRegion->mv[0]    = 0;	//////////////////////////////////
        pRegion->mv[1]    = 0;	//++ 设置每个8*8块的运动向量
        pRegion->mv[2]    = 0;	//////////////////////////////////
      }
      else
      {
        ii              = 4*mbx + (i%2)*2;// + BLOCK_SIZE;
        jj              = 4*mby + (i/2)*2;
        if (currMB->b8mode[i]>=5 && currMB->b8mode[i]<=7)  // SMALL BLOCKS	//++ 5：8x4，6：4x8，7：4x4
        {
          pRegion->mv[0]  = (dec_picture->mv[LIST_0][ii][jj][0] + dec_picture->mv[LIST_0][ii+1][jj][0] + dec_picture->mv[LIST_0][ii][jj+1][0] + dec_picture->mv[LIST_0][ii+1][jj+1][0] + 2)/4;
          pRegion->mv[1]  = (dec_picture->mv[LIST_0][ii][jj][1] + dec_picture->mv[LIST_0][ii+1][jj][1] + dec_picture->mv[LIST_0][ii][jj+1][1] + dec_picture->mv[LIST_0][ii+1][jj+1][1] + 2)/4;
        }
        else // 16x16, 16x8, 8x16, 8x8
        {
          pRegion->mv[0]  = dec_picture->mv[LIST_0][ii][jj][0];
          pRegion->mv[1]  = dec_picture->mv[LIST_0][ii][jj][1];
//          pRegion->mv[0]  = dec_picture->mv[LIST_0][4*mbx+(i%2)*2+BLOCK_SIZE][4*mby+(i/2)*2][0];
//          pRegion->mv[1]  = dec_picture->mv[LIST_0][4*mbx+(i%2)*2+BLOCK_SIZE][4*mby+(i/2)*2][1];
        }
        erc_mvperMB      += mabs(pRegion->mv[0]) + mabs(pRegion->mv[1]);
        pRegion->mv[2]    = dec_picture->ref_idx[LIST_0][ii][jj];
      }
    }
  }
  else  //B-frame
  {
    for (i=0; i<4; i++)
    {
      ii                  = 4*mbx + (i%2)*2;// + BLOCK_SIZE;
      jj                  = 4*mby + (i/2)*2;
      pRegion             = currRegion + i;
      pRegion->regionMode = (currMB->mb_type  ==I16MB  ? REGMODE_INTRA      :
                             currMB->b8mode[i]==IBLOCK ? REGMODE_INTRA_8x8  : REGMODE_INTER_PRED_8x8);
      if (currMB->mb_type==I16MB || currMB->b8mode[i]==IBLOCK)  // INTRA
      {
        pRegion->mv[0]    = 0;
        pRegion->mv[1]    = 0;
        pRegion->mv[2]    = 0;
      }
      else
      {
        int idx = (dec_picture->ref_idx[0][ii][jj]<0)?1:0;
//        int idx = (currMB->b8mode[i]==0 && currMB->b8pdir[i]==2 ? LIST_0 : currMB->b8pdir[i]==1 ? LIST_1 : LIST_0);
//        int idx = currMB->b8pdir[i]==0 ? LIST_0 : LIST_1;
        mv                = dec_picture->mv[idx];
        pRegion->mv[0]    = (mv[ii][jj][0] + mv[ii+1][jj][0] + mv[ii][jj+1][0] + mv[ii+1][jj+1][0] + 2)/4;
        pRegion->mv[1]    = (mv[ii][jj][1] + mv[ii+1][jj][1] + mv[ii][jj+1][1] + mv[ii+1][jj+1][1] + 2)/4;
        erc_mvperMB      += mabs(pRegion->mv[0]) + mabs(pRegion->mv[1]);

        pRegion->mv[2]  = (dec_picture->ref_idx[idx][ii][jj]);
/*        
        if (currMB->b8pdir[i]==0 || (currMB->b8pdir[i]==2 && currMB->b8mode[i]!=0)) // forward or bidirect
        {
          pRegion->mv[2]  = (dec_picture->ref_idx[LIST_0][ii][jj]);
          ///???? is it right, not only "img->fw_refFrArr[jj][ii-4]"
        }
        else
        {
          pRegion->mv[2]  = (dec_picture->ref_idx[LIST_1][ii][jj]);
//          pRegion->mv[2]  = 0;
        }
        */
      }
    }
  }
}

/*!
 ************************************************************************
 * \brief
 *    set defaults for old_slice
 *    NAL unit of a picture"
 ************************************************************************
 */
void init_old_slice()
{
  old_slice.field_pic_flag = 0;

  old_slice.pps_id = INT_MAX;

  old_slice.frame_num = INT_MAX;

  old_slice.nal_ref_idc = INT_MAX;
  
  old_slice.idr_flag = 0;

  old_slice.pic_oder_cnt_lsb          = UINT_MAX;
  old_slice.delta_pic_oder_cnt_bottom = INT_MAX;

  old_slice.delta_pic_order_cnt[0] = INT_MAX;
  old_slice.delta_pic_order_cnt[1] = INT_MAX;

}

/*!
 ************************************************************************
 * \brief
 *    save slice parameters that are needed for checking of "first VCL
 *    NAL unit of a picture"
 ************************************************************************
 */
void exit_slice()
{

  old_slice.pps_id = img->currentSlice->pic_parameter_set_id;

  old_slice.frame_num = img->frame_num;

  old_slice.field_pic_flag = img->field_pic_flag;

  if(img->field_pic_flag)
  {
    old_slice.bottom_field_flag = img->bottom_field_flag;
  }

  old_slice.nal_ref_idc   = img->nal_reference_idc;
  
  old_slice.idr_flag = img->idr_flag;
  if (img->idr_flag)
  {
    old_slice.idr_pic_id = img->idr_pic_id;
  }

  if (active_sps->pic_order_cnt_type == 0)
  {
    old_slice.pic_oder_cnt_lsb          = img->pic_order_cnt_lsb;
    old_slice.delta_pic_oder_cnt_bottom = img->delta_pic_order_cnt_bottom;
  }

  if (active_sps->pic_order_cnt_type == 1)
  {
    old_slice.delta_pic_order_cnt[0] = img->delta_pic_order_cnt[0];
    old_slice.delta_pic_order_cnt[1] = img->delta_pic_order_cnt[1];
  }
}

/*!
 ************************************************************************
 * \brief
 *    detect if current slice is "first VCL NAL unit of a picture"
 ************************************************************************
 */
int is_new_picture()
{
  int result=0;

  result |= (old_slice.pps_id != img->currentSlice->pic_parameter_set_id);

  result |= (old_slice.frame_num != img->frame_num);

  result |= (old_slice.field_pic_flag != img->field_pic_flag);

  if(img->field_pic_flag && old_slice.field_pic_flag)
  {
    result |= (old_slice.bottom_field_flag != img->bottom_field_flag);
  }

  result |= (old_slice.nal_ref_idc   != img->nal_reference_idc);
  
  result |= ( old_slice.idr_flag != img->idr_flag);

  if (img->idr_flag && old_slice.idr_flag)
  {
    result |= (old_slice.idr_pic_id != img->idr_pic_id);
  }

  if (active_sps->pic_order_cnt_type == 0)
  {
    result |=  (old_slice.pic_oder_cnt_lsb          != img->pic_order_cnt_lsb);
    result |=  (old_slice.delta_pic_oder_cnt_bottom != img->delta_pic_order_cnt_bottom);
  }

  if (active_sps->pic_order_cnt_type == 1)
  {
    result |= (old_slice.delta_pic_order_cnt[0] != img->delta_pic_order_cnt[0]);
    result |= (old_slice.delta_pic_order_cnt[1] != img->delta_pic_order_cnt[1]);
  }

  return result;
}


/*!
 ************************************************************************
 * \brief
 *    decodes one slice
 ************************************************************************
 */
/*decode_one_slice 函数才真正开始解码,看一下 read_one_macroblock 函数*/
void decode_one_slice(struct img_par *img,struct inp_par *inp) //片解码函数
{

  Boolean end_of_slice = FALSE;
  int read_flag;
  img->cod_counter=-1;

  set_ref_pic_num(); //首先用set_ref_pic_num()函数初始化参考帧序列，之后进入循环，对每一个宏块解码

  if (img->type == B_SLICE)
      compute_collocated(Co_located, listX);	//++ 为直接预测模式做一些准备工作：获取co_located图像、计算mv_scale

  //reset_ec_flags();

  while (end_of_slice == FALSE) // loop over macroblocks
  {

#if TRACE
  fprintf(p_trace,"\n*********** POC: %i (I/P) MB: %i Slice: %i Type %d **********\n", img->ThisPOC, img->current_mb_nr, img->current_slice_nr, img->type);
#endif

    // Initializes the current macroblock
    start_macroblock(img,inp, img->current_mb_nr);//初始化当前宏块
    // Get the syntax elements from the NAL
    read_flag = read_one_macroblock(img,inp);
    /*showcoeff(img);

	//++ 熵解码：包括解出宏块类型、预测模式、MVD、CBP(eodedBlock pattern)、残差（包括反量化操作）等

    /****************************************************************************************************************
	decode_one_macroblock()首先针对亮度块解码，之后对色度块解码。首先判断宏块是否采用16x16帧内预测，如果是，则通过函
	数intrapred_luma_16x16()解码。然后对一个16x16宏块中的每个4x4块分别解码（每个块的解码顺序见毕厚杰书的图6.49（P118）
	），直到所有块解码完成为止。解码4x4块时，首先判断是否为帧内编码，如果是，则通过函数intrapred对4x4块帧内解码。如果
	是帧间编码，则还要判断是P或是B帧。找到参考帧的匹配块后，用get_block()像素内插恢复像素值。itrans()对当前块的残差进
	行反变换，然后再将反变换值和预测值相加得到重构值。接着，将该重构值赋给dec_picture->imagY保存。其中de_picture是一
	个用来存放解码值的解构体，其成员imgy存放亮度解码值。最后是色度块解码。色度块解码和亮度块解码流程基本相同
	****************************************************************************************************************/
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//if(IS_OLDINTRA (&img->mb_data[img->current_mb_nr]))
	//  decryption_intra_mode(img);
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	decode_one_macroblock(img,inp);				//++ 反变换及运动补偿：反量化反变换、运动补偿、像素重构等
    if(img->MbaffFrameFlag && dec_picture->mb_field[img->current_mb_nr])
    {
      img->num_ref_idx_l0_active >>= 1;
      img->num_ref_idx_l1_active >>= 1;
    }

    ercWriteMBMODEandMV(img,inp);	//++ 写入各个8*8块的预测模式及运动向量到错误隐藏变量中

    end_of_slice=exit_macroblock(img,inp,(!img->MbaffFrameFlag||img->current_mb_nr%2));
  }

  exit_slice();
  //reset_ec_flags();

}


void decode_slice(struct img_par *img,struct inp_par *inp, int current_header)
//decode_slice 函数是片解码函数，对之前通过函数read_new_slice()读取的数据片进行解码,
//该函数中最重要的是 decode_one_slice 函数(这才真正开始解码)，再进入 
{
  Slice *currSlice = img->currentSlice;

  if (active_pps->entropy_coding_mode_flag)
  {
    init_contexts (img);
    cabac_new_slice();
  }

  if ( (active_pps->weighted_bipred_idc > 0  && (img->type == B_SLICE)) || (active_pps->weighted_pred_flag && img->type !=I_SLICE))
    fill_wp_params(img);

  //printf("frame picture %d %d %d\n",img->structure,img->ThisPOC,img->direct_type);
  

  // decode main slice information
  if ((current_header == SOP || current_header == SOS) && currSlice->ei_flag == 0)
    decode_one_slice(img,inp);//片解码函数
    
  // setMB-Nr in case this slice was lost
//  if(currSlice->ei_flag)  
//    img->current_mb_nr = currSlice->last_mb_nr + 1;

}


/*!
 ************************************************************************
 * \brief
 *    Prepare field and frame buffer after frame decoding
 ************************************************************************
 */
void frame_postprocessing(struct img_par *img, struct inp_par *inp)
{
}

/*!
 ************************************************************************
 * \brief
 *    Prepare field and frame buffer after field decoding
 ************************************************************************
 */
void field_postprocessing(struct img_par *img, struct inp_par *inp)
{
  img->number /= 2;
}



void reset_wp_params(struct img_par *img)
{
  int i,comp;
  int log_weight_denom;

  for (i=0; i<MAX_REFERENCE_PICTURES; i++)
  {
    for (comp=0; comp<3; comp++)
    {
      log_weight_denom = (comp == 0) ? img->luma_log2_weight_denom : img->chroma_log2_weight_denom;
      img->wp_weight[0][i][comp] = 1<<log_weight_denom;
      img->wp_weight[1][i][comp] = 1<<log_weight_denom;
    }
  }
}


void fill_wp_params(struct img_par *img)
{
  int i, j, k;
  int comp;
  int log_weight_denom;
  int p0, pt;
  int bframe = (img->type==B_SLICE);
  int max_bwd_ref, max_fwd_ref;
  int x,z;

  max_fwd_ref = img->num_ref_idx_l0_active;
  max_bwd_ref = img->num_ref_idx_l1_active;

  if (active_pps->weighted_bipred_idc == 2 && bframe)
  {
    img->luma_log2_weight_denom = 5;
    img->chroma_log2_weight_denom = 5;
    img->wp_round_luma = 16;
    img->wp_round_chroma = 16;

    for (i=0; i<MAX_REFERENCE_PICTURES; i++)
    {
      for (comp=0; comp<3; comp++)
      {
        log_weight_denom = (comp == 0) ? img->luma_log2_weight_denom : img->chroma_log2_weight_denom;
        img->wp_weight[0][i][comp] = 1<<log_weight_denom;
        img->wp_weight[1][i][comp] = 1<<log_weight_denom;
        img->wp_offset[0][i][comp] = 0;
        img->wp_offset[1][i][comp] = 0;
      }
    }
  }

  if (bframe)
  {
    for (i=0; i<max_fwd_ref; i++)
    {
      for (j=0; j<max_bwd_ref; j++)
      {
        for (comp = 0; comp<3; comp++)
        {
          log_weight_denom = (comp == 0) ? img->luma_log2_weight_denom : img->chroma_log2_weight_denom;
          if (active_pps->weighted_bipred_idc == 1)
          {
            img->wbp_weight[0][i][j][comp] =  img->wp_weight[0][i][comp];
            img->wbp_weight[1][i][j][comp] =  img->wp_weight[1][j][comp];
          }
          else if (active_pps->weighted_bipred_idc == 2)
          {
            pt = listX[LIST_1][j]->poc - listX[LIST_0][i]->poc;
            if (pt == 0 || listX[LIST_1][j]->is_long_term || listX[LIST_0][i]->is_long_term)
            {
              img->wbp_weight[0][i][j][comp] =   32;
              img->wbp_weight[1][i][j][comp] =   32;
            }
            else
            {
              p0 = img->ThisPOC - listX[LIST_0][i]->poc;
              
              x = (16384 + (pt>>1))/pt;
              z = Clip3(-1024, 1023, (x*p0 + 32 )>>6);
              
              img->wbp_weight[1][i][j][comp] = z >> 2;
              img->wbp_weight[0][i][j][comp] = 64 - img->wbp_weight[1][i][j][comp];
              if (img->wbp_weight[1][i][j][comp] < -64 || img->wbp_weight[1][i][j][comp] > 128)
              {
                img->wbp_weight[0][i][j][comp] = 32;
                img->wbp_weight[1][i][j][comp] = 32;
                img->wp_offset[0][i][comp] = 0;
                img->wp_offset[1][i][comp] = 0;
              }
            }
          }
        }
      }
    }
  }

  if (bframe && img->MbaffFrameFlag)
  {
    for (i=0; i<2*max_fwd_ref; i++)
    {
      for (j=0; j<2*max_bwd_ref; j++)
      {
        for (comp = 0; comp<3; comp++)
        {
          for (k=2; k<6; k+=2)
          {
            img->wp_offset[k+0][i][comp] = img->wp_offset[0][i/2][comp];
            img->wp_offset[k+1][i][comp] = img->wp_offset[1][i/2][comp];

            log_weight_denom = (comp == 0) ? img->luma_log2_weight_denom : img->chroma_log2_weight_denom;
            if (active_pps->weighted_bipred_idc == 1)
            {
              img->wbp_weight[k+0][i][j][comp] =  img->wp_weight[0][i/2][comp];
              img->wbp_weight[k+1][i][j][comp] =  img->wp_weight[1][j/2][comp];
            }
            else if (active_pps->weighted_bipred_idc == 2)
            {
              pt = listX[k+LIST_1][j]->poc - listX[k+LIST_0][i]->poc;
              if (pt == 0 || listX[k+LIST_1][j]->is_long_term || listX[k+LIST_0][i]->is_long_term)
              {
                img->wbp_weight[k+0][i][j][comp] =   32;
                img->wbp_weight[k+1][i][j][comp] =   32;
              }
              else
              {
                p0 = ((k==2)?img->toppoc:img->bottompoc) - listX[k+LIST_0][i]->poc;
                
                x = (16384 + (pt>>1))/pt;
                z = Clip3(-1024, 1023, (x*p0 + 32 )>>6);

                img->wbp_weight[k+1][i][j][comp] = z >> 2;
                img->wbp_weight[k+0][i][j][comp] = 64 - img->wbp_weight[k+1][i][j][comp];
                if (img->wbp_weight[k+1][i][j][comp] < -64 || img->wbp_weight[k+1][i][j][comp] > 128)
                {
                  img->wbp_weight[k+1][i][j][comp] = 32;
                  img->wbp_weight[k+0][i][j][comp] = 32;
                  img->wp_offset[k+0][i][comp] = 0;
                  img->wp_offset[k+1][i][comp] = 0;
                }
              }
            }
          }
        }
      }
    }
  }
}
