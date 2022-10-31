/* Unpack a LIME formated file 
   Each LIME record payload becomes a file */
/* F. Maresca 12/22/04 based on code by Balint Joo */

/* Usage...

   lime_unpack <lime_file>

   Files are unpacked into a directory named after the lime_file with
   names that encode the message number, record number, and LIME type.

   lime_file.contents/msgnn.recnn.lime_type

*/

#include <lime_config.h>
#include <stdio.h>
#include <lime.h>
#include <lime_fixed_types.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>

void Truncate( FILE *fp, LimeReader *reader, char *limefile, int Msg, off_t TruncateAt )
{
  limeDestroyReader(reader);
  fclose(fp);
  printf("Truncate after msg %d = %llu bytes\n", Msg, (unsigned long long)TruncateAt);
  truncate( limefile, TruncateAt );
  exit( EXIT_SUCCESS );
}

int main(int argc, char *argv[])
{
  FILE *fp;
  int TruncateAfter;
  char *limefile;
  char *lime_type;
  int rec, msg, status;
  off_t ThisRecord;
  LimeReader *reader;
  n_uint64_t nbytes;
  size_t bytes_pad;
  int MB_flag, ME_flag;
  
  if( argc < 3 )
    {
      fprintf(stderr, "Usage: %s <truncate_after> <lime_file>\n", argv[0]);
      return EXIT_FAILURE;
    }

  TruncateAfter = atoi( argv[1] );
  if( !TruncateAfter )
  {
    fprintf(stderr, "%s is not a record number\n", argv[1]);
    return EXIT_FAILURE;
  }
  limefile = argv[2];
  
  /* Open LIME file for reading */
  fp = DCAPL(fopen)(limefile, "r");
  if(fp == (FILE *)NULL) 
    {
      fprintf(stderr,"Unable to open file %s for reading\n", limefile);
      return EXIT_FAILURE;
    }
  
  /* Open LIME reader */
  reader = limeCreateReader(fp);
  if( reader == (LimeReader *)NULL ) 
    {
      fprintf(stderr, "Unable to open LimeReader\n");
      return EXIT_FAILURE;
    }
  
  /* Loop over LIME records */
  rec = 0;
  msg = 0;

  printf("   bytes file\n");

  ThisRecord = limeGetReaderPointer(reader);
  while( (status = limeReaderNextRecord(reader)) == LIME_SUCCESS )
  {
    rec++;
    
    nbytes    = limeReaderBytes(reader);
    lime_type = limeReaderType(reader);
    bytes_pad = limeReaderPadBytes(reader);
    MB_flag   = limeReaderMBFlag(reader);
    ME_flag   = limeReaderMEFlag(reader);
    
    if (MB_flag == 1)
    {
      rec = 1;
      ++msg;
      printf("\n");
    }

    if( msg - 1 == TruncateAfter )
    {
      Truncate( fp, reader, limefile, TruncateAfter, ThisRecord );
    }

    /* Announce file */
    printf("%d:%d %8llu %8llu %s\n", msg, rec,
           (unsigned long long)ThisRecord, (unsigned long long)nbytes, lime_type);
    ThisRecord = limeGetReaderPointer(reader);
  }
  if( status != LIME_EOF )
  {
    fprintf(stderr, "limeReaderNextRecord returned status = %d\n", status);
    return EXIT_FAILURE;
  }

  if( msg == TruncateAfter )
  {
    printf("File already contains %d messages\n", TruncateAfter);
    return EXIT_SUCCESS;
  }

  limeDestroyReader(reader);
  fclose(fp);

  fprintf(stderr, "EOF hit before message %d\n", TruncateAfter);
  return EXIT_FAILURE;
}

