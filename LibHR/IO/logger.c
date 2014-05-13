/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *   
* All rights reserved.                                                      * 
\***************************************************************************/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <stdarg.h>
#include "logger.h"
#include "global.h"

/* 
 * ***********************************************
 * Simple output logging facility
 * ***********************************************
 */

typedef struct _record {
	char *name;
	FILE *file;
	struct _record *next;
} record;

typedef struct _lrecord {
	char *name;
	int level;
	struct _lrecord *next;
} lrecord;

static record *filemap=0; /* list the mappings to files */
static record *files=0; /* list of open files */
static record *default_out=0; /* this is a special list with only one element if stdout is mapped */
static lrecord *levels=0; /* list of verbosity levels of IDs */
static int mapchanged=1; /* to keep track of changes to previous maps */

static int verblevel=0; /* the default verbosity level */

static void cleanup() {
	record *curr;
	lrecord *lcurr;
	void *tmp;

	/* reset default out to stdout */
	if(default_out!=0) {
		free(default_out->name);
		fclose(default_out->file);
		free(default_out);
		default_out=0;
	}

	/* free files list and close files */
	curr=files;
	while(curr!=0) {
		free(curr->name);
		fclose(curr->file);
		tmp=curr;
		curr=curr->next;
		free(tmp);
	}
	files=0;
	
	/* now free map list */
	curr=filemap;
	while(curr!=0) {
		free(curr->name);
		tmp=curr;
		curr=curr->next;
		free(tmp);
	}
	filemap=0;

	/* free the verbosity map */
	lcurr=levels;
	while(lcurr!=0) {
		free(lcurr->name);
		tmp=lcurr;
		lcurr=lcurr->next;
		free(tmp);
	}
	levels=0;

	mapchanged=1;
}

/* find a record with the given file in the list curr */
static record *findfile(record *curr, FILE *file) {

	assert(file!=0);
	while(curr!=0) {
		if(file==curr->file)
			return curr;
		curr=curr->next;
	}
	return 0;
}

/* find a record with the given name in the record list *curr */
static record *findname(record *curr, char *name) {

	assert(name!=0);
	while(curr!=0) {
		if(strcmp(name,curr->name)==0)
			return curr;
		curr=curr->next;
	}
	return 0;
}

/* find a record with the given name in the verbosity list */
static lrecord *lfindname(lrecord *curr, char *name) {

	assert(name!=0);
	while(curr!=0) {
		if(strcmp(name,curr->name)==0)
			return curr;
		curr=curr->next;
	}
	return 0;
}

/* create a new record and put it at the beginning of the given list */
static void addrecord(record **list, char* name, FILE* file) {
	record *new;

	assert(name!=0);
	assert(file!=0);

	new=malloc(sizeof(*new));
	new->name=malloc(sizeof(char)*(strlen(name)+1));
	strcpy(new->name,name);
	new->file=file;
	new->next=*list;
	*list=new;

	mapchanged=1;

}

/* create a new record and put it at the beginning of the given list */
static void addlrecord(lrecord **list, char* name, int level) {
	lrecord *new;

	assert(name!=0);

	new=malloc(sizeof(*new));
	new->name=malloc(sizeof(char)*(strlen(name)+1));
	strcpy(new->name,name);
	new->level=level;
	new->next=*list;
	*list=new;

	mapchanged=1;

}

/* remove a record from a list */
static void delrecord(record **list, record *rd) {
	record *cur=(*list)->next;
	record *pcur=*list;

	assert(rd!=0);

	mapchanged=1;

	free(rd->name);

	if (rd==*list) { /* rd id the first record */
		*list=(*list)->next;
		free(rd);
		return;
	}
	do {
		if(rd==cur){
			pcur->next=rd->next;
			free(rd);
			return;
		}
		pcur=cur;
		cur=cur->next;
	} while(cur!=0);
		
	assert(0); /* this means the the record wasn't in the list */

}

/* remove a record from a list */
static void dellrecord(lrecord **list, lrecord *rd) {
	lrecord *cur=(*list)->next;
	lrecord *pcur=*list;

	assert(rd!=0);

	mapchanged=1;

	free(rd->name);

	if (rd==*list) { /* rd id the first record */
		*list=(*list)->next;
		free(rd);
		return;
	}
	do {
		if(rd==cur){
			pcur->next=rd->next;
			free(rd);
			return;
		}
		pcur=cur;
		cur=cur->next;
	} while(cur!=0);
		
	assert(0); /* this means the the record wasn't in the list */

}

int logger_unmap(char *name){
	record *rd;
	FILE *file;

	if(name==0)
		return 1;

	rd=findname(filemap,name);
	if(rd==0) /* record do not exist: nothing to do */
		return 0;

	file=rd->file;
	delrecord(&filemap,rd);

	if(findfile(filemap,file)==0) {
		/* no other mapping to this file exists */
		/* close file */
		rd=findfile(files,file);
		assert(rd!=0);
		fclose(file);
		delrecord(&files,rd);
	}

	return 0;

}

/* reset mappping */
int logger_reset() {
	cleanup();
	verblevel=0;

	return 0;
}

/* link the ID name to the file with name filename */
/* Returns:
 * 0 => success
 * 1 => invalid name
 * 2 => invalid filename
 * 3 => cannot open new file
 */
int logger_map(char *name, char *filename) {
	static int init=1;
	char openmode='w';
	FILE *newfd;
	record *rd;

	if(name==0)
		return 1;
	if(filename==0)
		return 2;

	/* register cleanup function */
	if (init) {
		atexit(&cleanup);
		init=0;
	}

	/* check if we want append to file */
	/* the filename must begin with >> */
	/* if the filename start with > this is ignored */
	if(filename[0]=='>'){
		++filename;
		if(filename[0]=='>'){
			++filename;
			openmode='a';
		}
	}

	/* check if filename is the default_out of the logger */
	if(default_out!=0){
		if(strcmp(default_out->name,filename)==0){
			logger_unmap(name);
			return 0; /* mapping is not necessary */
		}
	}

	/* check if name is already mapped to filename */
	rd=findname(filemap,name);
	if (rd!=0) {
		rd=findfile(files,rd->file);
		if(strcmp(filename,rd->name)==0) {
			return 0;
		}
	}

	rd=findname(files,filename);
	if (rd==0) {
		/* open new file */
		newfd=fopen(filename,&openmode);
		if(newfd==0) {
			return 3;
		}
		addrecord(&files,filename,newfd);
	} else {
		newfd=rd->file;
	}

	logger_unmap(name);

	addrecord(&filemap,name,newfd);

	return 0;

}

void logger_setlevel(char *name, int v){
	lrecord *rd;

	if(name==0) { /* set the default verbosity level */
		verblevel=v;
		return;
	}

	rd=lfindname(levels,name);
	if(rd!=0){
		rd->level=v;
		return;
	}

	addlrecord(&levels,name,v);

}

void logger_set_input(input_logger *logger){
  if(logger->def_log_lvl==-1)
    logger->def_log_lvl=10;
  
  logger_setlevel(0,logger->def_log_lvl);


  if(logger->inverter_log_lvl==-1)
    logger->inverter_log_lvl=logger->def_log_lvl;
  else {
    logger_setlevel("INVERTER",logger->inverter_log_lvl);
  }
  if(logger->forcestat_log_lvl==-1)
    logger->forcestat_log_lvl=logger->def_log_lvl;
  else {
    logger_setlevel("FORCE-STAT",logger->forcestat_log_lvl);
  }
}

int logger_getlevel(char *name){
	lrecord *rd;
	if (name==0) { /* default verblevel */
		return verblevel;
	}

	rd=lfindname(levels,name);
	if(rd!=0){
		return rd->level;
	}

	/* if no record exist return default verbosity level */
	return verblevel;

}

/* this function reset the verbosity level of name to default */
void logger_rmlevel(char *name){
	lrecord *rd;

	if(name==0)
		return; /* nothing to do */

	rd=lfindname(levels,name);
	if(rd==0) 
		return; /* no record found: already at defaut level */

	dellrecord(&levels,rd);
	
}

int logger_stdout(char *filename) {
	FILE *fd;
	char openmode='w';
	record *rd;

	mapchanged=1; /* assume we are about to change the map */

	/* reset default_out to zero if filename=0 */
	if(filename==0){
		if(default_out!=0) {
			free(default_out->name);
			fclose(default_out->file);
			free(default_out);
			default_out=0;
		}
		/* reset default_out to zero */
		return 0;
	}
	
	/* if filename start with >> open the file in open mode */
	if(filename[0]=='>'){
		++filename;
		if(filename[0]=='>'){
			++filename;
			openmode='a';
		}
	}

	/* check if filename is already open */
  rd=findname(files,filename);
	if(rd!=0) {
		fd=rd->file;
		/* we must delete all entries in filemap with file=filename */
		delrecord(&files,rd);
		while((rd=findfile(filemap,fd))!=0)
			delrecord(&filemap,rd);
	} else {
		/* try to open the new file */
		fd=fopen(filename,&openmode);
		if(fd==0) {
			mapchanged=0; /* map has NOT be changed */
			return 1;
		}
	}

	/* change logger stdout */
	if(default_out!=0) {
		fclose(default_out->file);
		delrecord(&default_out,default_out);
	}

	addrecord(&default_out,filename,fd);

	return 0;
}

static void mycpyname(char **dst, char *src){
	*((*dst)++)='[';
	while(*src) {
		*((*dst)++)=*(src++);
	}
	*((*dst)++)=']';
}

static int mycpytonl(char **dst, char **src){
	while(**src) {
		*((*dst)++)=**src;
		if(*((*src)++)=='\n') {
			if(**src=='\0')
				return 0;
			else 
				return 1;
		}
	}
	return 0;
}

static int logger_inactive=0;

void logger_enable() {
  logger_inactive=0;
}

void logger_disable() {
  logger_inactive=1;
}

int lprintf(char *name, int level, char *format, ...) {
	va_list args;
	static record *lastrec=0;
	static char lastname[512]={0};
	static FILE *lastfd=0;
	static char buf[1024]; 
	static char alevel[16];
	char *cur=&buf[0];
	int ret;
	lrecord *vrd;
	static int lastvlevel;
	static int newline=1;
	int islast=1;

	if(logger_inactive || name==0) /* no name: print nothing and return and error */
		return -1;

	/* compare current name with last name if map has not changed */
	if(mapchanged || strcmp(name,lastname)!=0) {
		islast=0;
		mapchanged=0;
		lastrec=findname(filemap,name);
		if(lastrec==0){
			lastfd=(default_out==0)?stdout:default_out->file;
		} else {
			lastfd=lastrec->file;
		}
		strcpy(&lastname[0],name);
		vrd=lfindname(levels,name);
		lastvlevel=verblevel;
		if(vrd!=0)
			lastvlevel=vrd->level;
	}

	/* check verbosity level */
	if(lastvlevel<level)
		return 0;

	va_start(args, format);

	sprintf(alevel,"%d",level);
	if(newline) {
		mycpyname(&cur,name);
		mycpyname(&cur,alevel);
	} else if (!islast) {
		*(cur++)='\n';
		mycpyname(&cur,name);
		mycpyname(&cur,alevel);
	}
	while(mycpytonl(&cur,&format)){
		mycpyname(&cur,name);
		mycpyname(&cur,alevel);
	}
	newline=(*(cur-1)=='\n')?1:0;
	*cur='\0';

	ret=vfprintf(lastfd,&buf[0],args);
#ifdef IO_FLUSH
	fflush(lastfd);
#endif

	va_end(args);

	return ret;
	
}



