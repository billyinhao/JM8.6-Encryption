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
 *  \mainpage
 *     This is the H.264/AVC decoder reference software. For detailed documentation
 *     see the comments in each file.
 *
 *  \author
 *     The main contributors are listed in contributors.h
 *
 *  \version
 *     JM 8.6
 *
 *  \note
 *     tags are used for document system "doxygen"
 *     available at http://www.doxygen.org
 */
/*!
 *  \file
 *     ldecod.c
 *  \brief
 *     H.264/AVC reference decoder project main()
 *  \author
 *     Main contributors (see contributors.h for copyright, address and affiliation details)
 *     - Inge Lille-Lang鴜       <inge.lille-langoy@telenor.com>
 *     - Rickard Sjoberg         <rickard.sjoberg@era.ericsson.se>
 *     - Stephan Wenger          <stewe@cs.tu-berlin.de>
 *     - Jani Lainema            <jani.lainema@nokia.com>
 *     - Sebastian Purreiter     <sebastian.purreiter@mch.siemens.de>
 *     - Byeong-Moon Jeon        <jeonbm@lge.com>
 *     - Gabi Blaettermann       <blaetter@hhi.de>
 *     - Ye-Kui Wang             <wyk@ieee.org>
 *     - Karsten Suehring        <suehring@hhi.de>
 *
 ***********************************************************************
 */

#include "contributors.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <sys/timeb.h>

#if defined WIN32
  #include <conio.h>
#endif

#include <assert.h>

#include "global.h"
#include "rtp.h"
#include "memalloc.h"
#include "mbuffer.h"
#include "leaky_bucket.h"
#include "fmo.h"
#include "annexb.h"
#include "output.h"
#include "cabac.h"

#include "erc_api.h"

#define JM          "8"
#define VERSION     "8.6"

#define LOGFILE     "log.dec"
#define DATADECFILE "dataDec.txt"
#define TRACEFILE   "trace_dec.txt"

extern objectBuffer_t *erc_object_list;
extern ercVariables_t *erc_errorVar;
extern ColocatedParams *Co_located;

// I have started to move the inp and img structures into global variables.
// They are declared in the following lines.  Since inp is defined in conio.h
// and cannot be overridden globally, it is defined here as input
//
// Everywhere, input-> and img-> can now be used either globally or with
// the local override through the formal parameter mechanism

extern FILE* bits;
extern StorablePicture* dec_picture;

struct inp_par    *input;       //!< input parameters from input configuration file
struct snr_par    *snr;         //!< statistics
struct img_par    *img;         //!< image parameters
//////////////////////////////////////////////////////////////
encryption_parameter en_ps;
/////////////////////////////////////////////////////////////

int global_init_done = 0;
/////////////////////////////////////////////////////////////////////////////////
double u=4, init_con=0.325;
int cn=8;
void pre_process() {
     
    // 基于 Logistic 映射，通过迭代产生模拟混沌序列
    // 热身运算 200 次，使其完全进入混沌状态
    int i; 
    for(i=0; i<200; i++) {
        init_con = u*init_con*(1-init_con);
        if(init_con == 0.5) // 避免产生后续全 1 的序列 
            init_con = 0.4;
    }        
}
void output_encryption_information(encryption_parameter * en_p)
{
  char c;
  short inte;
  printf("Do you want to activate the fuction of decryption?(y/n/s):");
  c=getchar();
  while(c!='y'&&c!='n'&&c!='s')
  {
   printf("Can't identify %c,please try again.\n",c);
   while(getchar()!='\n');
   c=getchar();
  }
   if(c=='n')
   {
    en_p->intra_mode=0;
    en_p->chroma_mode=0;
    en_p->DCT=0;
    en_p->MVD=0;
   }
   if(c=='y')
   {
    printf("Please enter the intensity of encryption you want to get(1-6):");
    scanf("%hd",&inte);
    while(inte<1||inte>6)
    {
     printf("%hd scale out the range(1-6),please try again.\n",inte);
     scanf("%hd",&inte);
    }
    switch(inte)
    {
     case 1:
     en_p->intra_mode=0;
     en_p->chroma_mode=0;
     en_p->DCT=1;
     en_p->MVD=0;
     break;
     case 2:
     en_p->intra_mode=1;
     en_p->chroma_mode=1;
     en_p->DCT=0;
     en_p->MVD=1;
     break;
     case 3:
     en_p->intra_mode=1;
     en_p->chroma_mode=1;
     en_p->DCT=1;
     en_p->MVD=1;
     break;
     case 4:
     en_p->intra_mode=1;
     en_p->chroma_mode=1;
     en_p->DCT=2;
     en_p->MVD=1;
     break;
     case 5:
     en_p->intra_mode=1;
     en_p->chroma_mode=1;
     en_p->DCT=3;
     en_p->MVD=1;
     break;
     case 6:
     en_p->intra_mode=1;
     en_p->chroma_mode=1;
     en_p->DCT=4;
     en_p->MVD=1;
     break;
    }
   }
   if(c=='s')
   {
    printf("Now set the encryption parameter\n");
    printf("intra_mode(0-2):");
    scanf("%d",&en_p->intra_mode);
    while(en_p->intra_mode<0||en_p->intra_mode>2)
    {
     printf("%d scale out the range(0-2),please try again.\n",en_p->intra_mode);
     scanf("%d",&en_p->intra_mode);
    }
     printf("chroma_mode(0-1):");
     scanf("%d",&en_p->chroma_mode);
    while(en_p->chroma_mode<0||en_p->chroma_mode>1)
    {
     printf("%d scale out the range(0-1),please try again.\n",en_p->chroma_mode);
     scanf("%d",&en_p->chroma_mode);
    }
     printf("DCT(0-4):");
     scanf("%d",&en_p->DCT);
    while(en_p->chroma_mode<0||en_p->chroma_mode>4)
    {
     printf("%d scale out the range(0-1),please try again.\n",en_p->DCT);
     scanf("%d",&en_p->DCT);
    }
     printf("MVD(0-2):");
     scanf("%d",&en_p->MVD);
    while(en_p->MVD<0||en_p->MVD>2)
    {
     printf("%d scale out the range(0-1),please try again.\n",en_p->MVD);
     scanf("%d",&en_p->MVD);
    }
   }
}
/////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////


/*!
 ***********************************************************************
 * \brief
 *    main function for TML decoder
 ***********************************************************************
 */

/******************************************************************************************************
首先应该清楚码流结构，H.264 码流是由一个个 NALU 组成的 ，而 NALU 又是由 NALU 头和 RBSP 数据组成 ，
这体现在标准文档中就是 7.3.1 小节 ，在第 7 章的 Syntax 表中，我们可以看见很多加粗了的字体，凡是加
粗了的字体都是需要编码的语法元素，7.3.1 第一个加粗的字体是 forbidden_zero_bit ，最后一列是 f(1)，
这说明该语法元素采用 1 比特定长编码，那么 forbidden_zero_bit 到底是什么意思呢 ？7.4.1 小节中做了解
释。7.3.1 第二个加粗的字体是 nal_ref_idc，该行最后一列是 u(2)，这说明该语法元素采用 u(2) 方式编码。
那么 u(2) 是什么编码方式呢？7.3 节上面的段落做了说明，对包括 ue、u、te等等编码方式都在那里做了说明。
7.3.1 第三个加粗的字体是 nal_unit_type，该语法元素编码方式跟前一个语法元素nal_ref_idc 是同样道理。
所有加粗的字体到底表示什么含义在 Semantics 一部分中都有相应说明。例如 7.3.1 小节中加粗字体的说明 就
是 7.4.1 小节中，7.3.2.1 小节中加粗的字体的说明就在 7.4.2.1 小节中，7.3 节与 7.4 节是完全对应的，7.3.1 
小节说明了一个 NALU 的结构，对于一个 H.264 码流中的一个 NALU ，该 NALU 的第一个比特必然是表示 
forbidden_zero_bit，该 NALU 的第二、三个比特必然是表示 nal_ref_idc，该 NALU 的第四~八个比特必然是表示
nal_unit_type，那么对于解码器，解码 H.264 码流必然要按顺序处理这些语法元素。看JM 解码器是如何按顺序处
理这些语法元素的，首先进入main函数的decode_one_frame 函数(decode_one_frame 以前的部分可以不做分析，因
为这些部分都是一些初始化工作，以及读文件操作，并没有解码)，再进入 read_new_slice 函数，read_new_slice 
函数的主要作用是读取 NALU 单元数据 ，附带做一些解码的准备工作 。在 read_new_slice 函数中我们最先遇到的
两个函数是 GetAnnexbNALU 和 GetRTPNALU，这就涉及到 H.264 的两种码流格式：GetAnnexbNALU 处理字节流格式的
码流  GetRTPNALU 处理 RTP 格式码流 
 *******************************************************************************************************/
 

int main(int argc, char **argv)
{
//////////////////////////////////////////////////////////////////////////
	 pre_process();
//////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////
	 output_encryption_information(&en_ps);
 ///////////////////////////////////////////////////////////////////////////////////////////////////////////

    // allocate memory for the structures
  if ((input =  (struct inp_par *)calloc(1, sizeof(struct inp_par)))==NULL) no_mem_exit("main: input");//被我注释了
  if ((snr =  (struct snr_par *)calloc(1, sizeof(struct snr_par)))==NULL) no_mem_exit("main: snr");//被我注释了
if ((img =  (struct img_par *)calloc(1, sizeof(struct img_par)))==NULL) no_mem_exit("main: img");//被我注释了

  // Read Configuration File
  if (argc != 2)
  {
    snprintf(errortext, ET_SIZE, "Usage: %s <config.dat> \n\t<config.dat> defines decoder parameters",argv[0]);//被我注释了
   error(errortext, 300);//被我注释了
  }

  init_conf(input, argv[1]);

  init_old_slice();
  
  switch (input->FileFormat)//input->FileFormat结构中存放着用户对文档格式设置的参数
  {
  case 0:
    OpenBitstreamFile (input->infile);//为“0”，表示比特流的格式.input->infile表示需要打开的视频文件
    break;
  case 1:
    OpenRTPFile (input->infile);//为“1”，表示RTP数据包格式.input->infile表示需要打开的视频文件
    break;
  default:
    printf ("Unsupported file format %d, exit\n", input->FileFormat);
  }
  // Allocate Slice data struct
  malloc_slice(input,img);

  init(img);

  dec_picture = NULL;

  dpb.init_done = 0;
  g_nFrame = 0;

//  init_dpb(input);
  init_out_buffer();

  img->idr_psnr_number=input->ref_offset;
  img->psnr_number=0;

  img->number=0;
  img->type = I_SLICE;
  img->dec_ref_pic_marking_buffer = NULL;

  // B pictures
  Bframe_ctr=0;

  // time for total decoding session
  tot_time = 0;
  while (decode_one_frame(img, input, snr) != EOS);//decode_one_frame 以前的部分都是一些初始化工作，
	                                               //以及读文件操作，并没有解码,这里才开始解码 
    

  report(input, img, snr);
  free_slice(input,img);
  FmoFinit();
  free_global_buffers();

  flush_dpb();

#ifdef PAIR_FIELDS_IN_OUTPUT
  flush_pending_output(p_out);
#endif

  CloseBitstreamFile();

  fclose(p_out);
//  fclose(p_out2);
  if (p_ref)
    fclose(p_ref);
#if TRACE
  fclose(p_trace);
#endif

  ercClose(erc_errorVar);

  free_dpb();
  uninit_out_buffer();

  free_collocated(Co_located);
  free (input);
  free (snr);
  free (img);
  
  //while( !kbhit() ); 
  return 0;
}


/*!
 ***********************************************************************
 * \brief
 *    Initilize some arrays
 ***********************************************************************
 */
void init(struct img_par *img)  //!< image parameters
{
  int i;

  // initilize quad matrix used in snr routine
  for (i=0; i <  256; i++)
  {
    img->quad[i]=i*i; // fix from TML 1, truncation removed
  }

  img->oldFrameSizeInMbs = -1;
}


/*!
 ************************************************************************
 * \brief
 *    Read input from configuration file
 *
 * \par Input:
 *    Name of configuration filename
 *
 * \par Output
 *    none
 ************************************************************************
 */
void init_conf(struct inp_par *inp,
               char *config_filename)
{
  FILE *fd;
  int NAL_mode;

  // read the decoder configuration file
  if((fd=fopen(config_filename,"r")) == NULL)
  {
    snprintf(errortext, ET_SIZE, "Error: Control file %s not found\n",config_filename);
    error(errortext, 300);
  }

  fscanf(fd,"%s",inp->infile);                // H.264 compressed input bitsream
  fscanf(fd,"%*[^\n]");

  fscanf(fd,"%s",inp->outfile);               // YUV 4:2:2 input format
  fscanf(fd,"%*[^\n]");

  fscanf(fd,"%s",inp->reffile);               // reference file
  fscanf(fd,"%*[^\n]");


  // Frame buffer size
  fscanf(fd,"%d,",&inp->dpb_size);   // may be overwritten in case of RTP NAL
  fscanf(fd,"%*[^\n]");
  if (inp->dpb_size < 1)
  {
    snprintf(errortext, ET_SIZE, "Decoded Picture Buffer Size is %d. It has to be at least 1",inp->dpb_size);
    error(errortext,1);
  }

  fscanf(fd,"%d",&(NAL_mode));                // NAL mode
    fscanf(fd,"%*[^\n]");

  switch(NAL_mode)
  {
  case 0:
    inp->FileFormat = PAR_OF_ANNEXB;
    break;
  case 1:
    inp->FileFormat = PAR_OF_RTP;
    break;
  default:
    snprintf(errortext, ET_SIZE, "NAL mode %i is not supported", NAL_mode);
    error(errortext,400);
  }

  fscanf(fd,"%d,",&inp->ref_offset);   // offset used for SNR computation
  fscanf(fd,"%*[^\n]");

  fscanf(fd,"%d,",&inp->poc_scale);   // offset used for SNR computation
  fscanf(fd,"%*[^\n]");


  if (inp->poc_scale < 1 || inp->poc_scale > 2)
  {
    snprintf(errortext, ET_SIZE, "Poc Scale is %d. It has to be 1 or 2",inp->poc_scale);
    error(errortext,1);
  }

#ifdef _LEAKYBUCKET_
  fscanf(fd,"%ld,",&inp->R_decoder);             // Decoder rate
  fscanf(fd, "%*[^\n]");
  fscanf(fd,"%ld,",&inp->B_decoder);             // Decoder buffer size
  fscanf(fd, "%*[^\n]");
  fscanf(fd,"%ld,",&inp->F_decoder);             // Decoder initial delay
  fscanf(fd, "%*[^\n]"); 
  fscanf(fd,"%s",inp->LeakyBucketParamFile);    // file where Leaky Bucket params (computed by encoder) are stored
  fscanf(fd,"%*[^\n]");
#endif

  fclose (fd);


#if TRACE
  if ((p_trace=fopen(TRACEFILE,"w"))==0)             // append new statistic at the end
  {
    snprintf(errortext, ET_SIZE, "Error open file %s!",TRACEFILE);
    error(errortext,500);
  }
#endif


  if ((p_out=fopen(inp->outfile,"wb"))==0)
  {
    snprintf(errortext, ET_SIZE, "Error open file %s ",inp->outfile);
    error(errortext,500);
  }
/*  if ((p_out2=fopen("out.yuv","wb"))==0)
  {
    snprintf(errortext, ET_SIZE, "Error open file %s ",inp->outfile);
    error(errortext,500);
  }*/

  fprintf(stdout,"--------------------------------------------------------------------------\n");
  fprintf(stdout," Decoder config file                    : %s \n",config_filename);
  fprintf(stdout,"--------------------------------------------------------------------------\n");
  fprintf(stdout," Input H.264 bitstream                  : %s \n",inp->infile);
  fprintf(stdout," Output decoded YUV 4:2:0               : %s \n",inp->outfile);
  fprintf(stdout," Output status file                     : %s \n",LOGFILE);
  if ((p_ref=fopen(inp->reffile,"rb"))==0)
  {
    fprintf(stdout," Input reference file                   : %s does not exist \n",inp->reffile);
    fprintf(stdout,"                                          SNR values are not available\n");
  }
  else
    fprintf(stdout," Input reference file                   : %s \n",inp->reffile);

  fprintf(stdout,"--------------------------------------------------------------------------\n");
#ifdef _LEAKYBUCKET_
  fprintf(stdout," Rate_decoder        : %8ld \n",inp->R_decoder);
  fprintf(stdout," B_decoder           : %8ld \n",inp->B_decoder);
  fprintf(stdout," F_decoder           : %8ld \n",inp->F_decoder);
  fprintf(stdout," LeakyBucketParamFile: %s \n",inp->LeakyBucketParamFile); // Leaky Bucket Param file
  calc_buffer(inp);
  fprintf(stdout,"--------------------------------------------------------------------------\n");
#endif
  fprintf(stdout,"POC must = frame# or field# for SNRs to be correct\n");
  fprintf(stdout,"Frame    POC   QP  SnrY    SnrU    SnrV   Time(ms)\n");
}

/*!
 ************************************************************************
 * \brief
 *    Reports the gathered information to appropriate outputs
 *
 * \par Input:
 *    struct inp_par *inp,
 *    struct img_par *img,
 *    struct snr_par *stat
 *
 * \par Output:
 *    None
 ************************************************************************
 */
void report(struct inp_par *inp, struct img_par *img, struct snr_par *snr)
{
  #define OUTSTRING_SIZE 255
  char string[OUTSTRING_SIZE];
  FILE *p_log;

#ifndef WIN32
  time_t  now;
  struct tm *l_time;
#else
  char timebuf[128];
#endif

  fprintf(stdout,"-------------------- Average SNR all frames ------------------------------\n");
  fprintf(stdout," SNR Y(dB)           : %5.2f\n",snr->snr_ya);
  fprintf(stdout," SNR U(dB)           : %5.2f\n",snr->snr_ua);
  fprintf(stdout," SNR V(dB)           : %5.2f\n",snr->snr_va);
  fprintf(stdout," Total decoding time : %.3f sec \n",tot_time*0.001);
  fprintf(stdout,"--------------------------------------------------------------------------\n");
  fprintf(stdout," Exit JM %s decoder, ver %s ",JM,VERSION);
  fprintf(stdout,"\n");
  // write to log file

  snprintf(string, OUTSTRING_SIZE, "%s", LOGFILE);
  if ((p_log=fopen(string,"r"))==0)                    // check if file exist
  {
    if ((p_log=fopen(string,"a"))==0)
    {
      snprintf(errortext, ET_SIZE, "Error open file %s for appending",string);
      error(errortext, 500);
    }
    else                                              // Create header to new file
    {
      fprintf(p_log," ------------------------------------------------------------------------------------------\n");
      fprintf(p_log,"|  Decoder statistics. This file is made first time, later runs are appended               |\n");
      fprintf(p_log," ------------------------------------------------------------------------------------------ \n");
      fprintf(p_log,"| Date  | Time  |    Sequence        |#Img|Format|SNRY 1|SNRU 1|SNRV 1|SNRY N|SNRU N|SNRV N|\n");
      fprintf(p_log," ------------------------------------------------------------------------------------------\n");
    }
  }
  else
  { 
    fclose(p_log);
    p_log=fopen(string,"a");                    // File exist,just open for appending
  }

#ifdef WIN32
  _strdate( timebuf );
  fprintf(p_log,"| %1.5s |",timebuf );

  _strtime( timebuf);
  fprintf(p_log," % 1.5s |",timebuf);
#else
  now = time ((time_t *) NULL); // Get the system time and put it into 'now' as 'calender time'
  time (&now);
  l_time = localtime (&now);
  strftime (string, sizeof string, "%d-%b-%Y", l_time);
  fprintf(p_log,"| %1.5s |",string );

  strftime (string, sizeof string, "%H:%M:%S", l_time);
  fprintf(p_log,"| %1.5s |",string );
#endif

  fprintf(p_log,"%20.20s|",inp->infile);

  fprintf(p_log,"%3d |",img->number);

  fprintf(p_log,"%6.3f|",snr->snr_y1);
  fprintf(p_log,"%6.3f|",snr->snr_u1);
  fprintf(p_log,"%6.3f|",snr->snr_v1);
  fprintf(p_log,"%6.3f|",snr->snr_ya);
  fprintf(p_log,"%6.3f|",snr->snr_ua);
  fprintf(p_log,"%6.3f|\n",snr->snr_va);

  fclose(p_log);

  snprintf(string, OUTSTRING_SIZE,"%s", DATADECFILE);
  p_log=fopen(string,"a");

  if(Bframe_ctr != 0) // B picture used
  {
    fprintf(p_log, "%3d %2d %2d %2.2f %2.2f %2.2f %5d "
      "%2.2f %2.2f %2.2f %5d "
      "%2.2f %2.2f %2.2f %5d %.3f\n",
      img->number, 0, img->qp,
      snr->snr_y1,
      snr->snr_u1,
      snr->snr_v1,
      0,
      0.0,
      0.0,
      0.0,
      0,
      snr->snr_ya,
      snr->snr_ua,
      snr->snr_va,
      0,
      (double)0.001*tot_time/(img->number+Bframe_ctr-1));
  }
  else
  {
    fprintf(p_log, "%3d %2d %2d %2.2f %2.2f %2.2f %5d "
      "%2.2f %2.2f %2.2f %5d "
      "%2.2f %2.2f %2.2f %5d %.3f\n",
      img->number, 0, img->qp,
      snr->snr_y1,
      snr->snr_u1,
      snr->snr_v1,
      0,
      0.0,
      0.0,
      0.0,
      0,
      snr->snr_ya,
      snr->snr_ua,
      snr->snr_va,
      0,
      (double)0.001*tot_time/img->number);
  }
  fclose(p_log);
}

/*!
 ************************************************************************
 * \brief
 *    Allocates a stand-alone partition structure.  Structure should
 *    be freed by FreePartition();
 *    data structures
 *
 * \par Input:
 *    n: number of partitions in the array
 * \par return
 *    pointer to DataPartition Structure, zero-initialized
 ************************************************************************
 */

DataPartition *AllocPartition(int n)
{
  DataPartition *partArr, *dataPart;
  int i;

  partArr = (DataPartition *) calloc(n, sizeof(DataPartition));
  if (partArr == NULL)
  {
    snprintf(errortext, ET_SIZE, "AllocPartition: Memory allocation for Data Partition failed");
    error(errortext, 100);
  }

  for (i=0; i<n; i++) // loop over all data partitions
  {
    dataPart = &(partArr[i]);
    dataPart->bitstream = (Bitstream *) calloc(1, sizeof(Bitstream));
    if (dataPart->bitstream == NULL)
    {
      snprintf(errortext, ET_SIZE, "AllocPartition: Memory allocation for Bitstream failed");
      error(errortext, 100);
    }
    dataPart->bitstream->streamBuffer = (byte *) calloc(MAX_CODED_FRAME_SIZE, sizeof(byte));	//++ MAX_CODED_FRAME_SIZE=8000000，一帧的字节数
    if (dataPart->bitstream->streamBuffer == NULL)
    {
      snprintf(errortext, ET_SIZE, "AllocPartition: Memory allocation for streamBuffer failed");
      error(errortext, 100);
    }
  }
  return partArr;
}




/*!
 ************************************************************************
 * \brief
 *    Frees a partition structure (array).  
 *
 * \par Input:
 *    Partition to be freed, size of partition Array (Number of Partitions)
 *
 * \par return
 *    None
 *
 * \note
 *    n must be the same as for the corresponding call of AllocPartition
 ************************************************************************
 */


void FreePartition (DataPartition *dp, int n)
{
  int i;

  assert (dp != NULL);
  assert (dp->bitstream != NULL);
  assert (dp->bitstream->streamBuffer != NULL);
  for (i=0; i<n; i++)
  {
    free (dp[i].bitstream->streamBuffer);
    free (dp[i].bitstream);
  }
  free (dp);
}


/*!
 ************************************************************************
 * \brief
 *    Allocates the slice structure along with its dependent
 *    data structures
 *
 * \par Input:
 *    Input Parameters struct inp_par *inp,  struct img_par *img
 ************************************************************************
 */
void malloc_slice(struct inp_par *inp, struct img_par *img)
{
  Slice *currSlice;

  img->currentSlice = (Slice *) calloc(1, sizeof(Slice));
  if ( (currSlice = img->currentSlice) == NULL)
  {
    snprintf(errortext, ET_SIZE, "Memory allocation for Slice datastruct in NAL-mode %d failed", inp->FileFormat);
    error(errortext,100);
  }
//  img->currentSlice->rmpni_buffer=NULL;
  //! you don't know whether we do CABAC hre, hence initialize CABAC anyway
  // if (inp->symbol_mode == CABAC)
  if (1)
  {
    // create all context models
    currSlice->mot_ctx = create_contexts_MotionInfo();
    currSlice->tex_ctx = create_contexts_TextureInfo();
  }
  currSlice->max_part_nr = 3;  //! assume data partitioning (worst case) for the following mallocs()
  currSlice->partArr = AllocPartition(currSlice->max_part_nr);
}


/*!
 ************************************************************************
 * \brief
 *    Memory frees of the Slice structure and of its dependent
 *    data structures
 *
 * \par Input:
 *    Input Parameters struct inp_par *inp,  struct img_par *img
 ************************************************************************
 */
void free_slice(struct inp_par *inp, struct img_par *img)
{
  Slice *currSlice = img->currentSlice;

  FreePartition (currSlice->partArr, 3);
//      if (inp->symbol_mode == CABAC)
  if (1)
  {
    // delete all context models
    delete_contexts_MotionInfo(currSlice->mot_ctx);
    delete_contexts_TextureInfo(currSlice->tex_ctx);
  }
  free(img->currentSlice);

  currSlice = NULL;
}

/*!
 ************************************************************************
 * \brief
 *    Dynamic memory allocation of frame size related global buffers
 *    buffers are defined in global.h, allocated memory must be freed in
 *    void free_global_buffers()
 *
 *  \par Input:
 *    Input Parameters struct inp_par *inp, Image Parameters struct img_par *img
 *
 *  \par Output:
 *     Number of allocated bytes
 ***********************************************************************
 */
int init_global_buffers()
{
  int memory_size=0;

  if (global_init_done)
  {
    free_global_buffers();
  }

  // allocate memory for reference frame in find_snr
  memory_size += get_mem2D(&imgY_ref, img->height, img->width);
  memory_size += get_mem3D(&imgUV_ref, 2, img->height_cr, img->width_cr);

  // allocate memory in structure img
  if(((img->mb_data) = (Macroblock *) calloc(img->FrameSizeInMbs, sizeof(Macroblock))) == NULL)
    no_mem_exit("init_global_buffers: img->mb_data");

  if(((img->intra_block) = (int*)calloc(img->FrameSizeInMbs, sizeof(int))) == NULL)
    no_mem_exit("init_global_buffers: img->intra_block");

  memory_size += get_mem2Dint(&(img->ipredmode), 4*img->PicWidthInMbs , 4*img->FrameHeightInMbs);

  memory_size += get_mem2Dint(&(img->field_anchor),4*img->FrameHeightInMbs, 4*img->PicWidthInMbs);

  memory_size += get_mem3Dint(&(img->wp_weight), 2, MAX_REFERENCE_PICTURES, 3);
  memory_size += get_mem3Dint(&(img->wp_offset), 6, MAX_REFERENCE_PICTURES, 3);
  memory_size += get_mem4Dint(&(img->wbp_weight), 6, MAX_REFERENCE_PICTURES, MAX_REFERENCE_PICTURES, 3);

  // CAVLC mem
  memory_size += get_mem3Dint(&(img->nz_coeff), img->FrameSizeInMbs, 4, 6);

  memory_size += get_mem2Dint(&(img->siblock),img->PicWidthInMbs  , img->FrameHeightInMbs);

  global_init_done = 1;

  img->oldFrameSizeInMbs = img->FrameSizeInMbs;

  return (memory_size);
}

/*!
 ************************************************************************
 * \brief
 *    Free allocated memory of frame size related global buffers
 *    buffers are defined in global.h, allocated memory is allocated in
 *    int init_global_buffers()
 *
 * \par Input:
 *    Input Parameters struct inp_par *inp, Image Parameters struct img_par *img
 *
 * \par Output:
 *    none
 *
 ************************************************************************
 */
void free_global_buffers()
{
  free_mem2D (imgY_ref);
  free_mem3D (imgUV_ref,2);

  // CAVLC free mem
  free_mem3Dint(img->nz_coeff, img->oldFrameSizeInMbs);

  free_mem2Dint(img->siblock);

  // free mem, allocated for structure img
  if (img->mb_data       != NULL) free(img->mb_data);

  free (img->intra_block);

  free_mem2Dint (img->ipredmode);

  free_mem2Dint(img->field_anchor);

  free_mem3Dint(img->wp_weight, 2);
  free_mem3Dint(img->wp_offset, 6);
  free_mem4Dint(img->wbp_weight, 6, MAX_REFERENCE_PICTURES);

  global_init_done = 0;

}



