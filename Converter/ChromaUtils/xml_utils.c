#include <string.h>
#include <stdlib.h>
#include <stdio.h>


void read_tag(char* from, char* to, char * tag)
{
  int ltag=strlen(tag)+1;
  char * addr_start=strstr(from,tag);
  
  if(addr_start ==NULL ) 
    {
      printf("Error no such tag(%s) in metadata\n",tag);
      exit(1);
    }
  addr_start+=ltag;


  char endtag[256]="/";
  strcat(endtag,tag);


  char * addr_end=strstr(from,endtag);
  
  if(addr_end ==NULL ) 
    {
      printf("Error no such tag(%s) in metadata\n",endtag);
      exit(1);
    }
  addr_end-=1;

  
  if(addr_end-addr_start>0)
    memcpy(to,addr_start,addr_end-addr_start);
  else
    {
      printf("Error badly formed metadata for tag(%s)\n",tag);
      exit(1);
    }
  
  
  to[addr_end-addr_start+1]='\0';
}

