#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "tifdata.h"
#include "tocgns.h"
#include "petsc.h"

extern FILE* yyin;
extern int yyparse();
extern int yylex();


int main(int argc, char ** args)
{
  char filename[64];
  PetscInitialize(&argc,&args,PETSC_NULL,PETSC_NULL);
  PetscTruth     flg;
  printf("------------------------------------------------------------\n");
  printf("TIF dump tool Ver 1.4 \n");
  printf("Extract mesh, solution or material information from TIF file.\n");
  printf("------------------------------------------------------------\n\n");

  PetscOptionsHasName(PETSC_NULL,"--help",&flg);

  if(argc<3 || flg)
  {
    printf("Help information:\n");
    printf("There are two kinds of TIF file. One for mesh and field data and\n");
    printf("another for 1D array data.\n\n");
    printf("Use:\n");
    printf(" dumptif tif_file -i                output TIF information\n");
    printf("to check TIF file information.\n");
    printf("\n");
    printf("For field data TIF file\n");
    printf(" dumptif tif_file -c cgns_file      save mesh and doping to CGNS file\n");
    printf(" dumptif tif_file -s solution_file  save field solution data to solution_file\n");
    printf(" dumptif tif_file -m material_file  save material information to material_file\n");
    printf("\n");
    printf("For 1D array data TIF file\n");
    printf(" dumptif tif_file -t                save 1D array data to file\n");
    printf("\n");
    exit(0);
  }

  yyin=fopen(args[1],"r");
  if(yyin==NULL)
  {
    printf("I can't open TIF file by filename %s\n",args[1]);
    exit(0);
  }
  if(yyparse()) exit(0);
  fclose(yyin);

  PetscOptionsHasName(PETSC_NULL,"-i",&flg);
  if(flg)
  {
    if(solhead.sol_num>0 && region_array.size()>0 && parameter_array.size()>0)
      printf("It is a field data TIF file. -c -s and -m option works.\n");
    if(transol.data_name_array.size()>0)
      printf("It is a 1D array data TIF file. Only -t option works.\n");
  }


  PetscOptionsGetString(PETSC_NULL,"-s",filename,63,&flg);
  if(flg)
  {
    if(!strlen(filename))
    {
      printf("You should give file name after -s option\n");
      exit(0);
    }
    if(solhead.sol_num>0)
    {
      FILE *fp = fopen(filename,"w");
      printf("There are %d solution data arrays in this TIF file,listed below\n",solhead.sol_num);
      for(int i=0;i<solhead.sol_num;i++)
        printf("\\%s  ",solhead.sol_name_array[i].c_str());

      printf("\n\nData will be written to %s as : x y soldata1 soldata2...",filename);

      for(size_t i=0;i<node_array.size();i++)
      {
        fprintf(fp,"%e\t%e\t",node_array[soldata[i].index].x,node_array[soldata[i].index].y);
        for(int j=0;j<solhead.sol_num;j++)
          fprintf(fp,"%e\t",soldata[i].data_array[j]);
        fprintf(fp,"\n");
      }
      fclose(fp);

      printf("finished\n");
    }
    else
    {
      printf("no solution data in this TIF file.\n");
    }
  }


  PetscOptionsGetString(PETSC_NULL,"-c",filename,63,&flg);
  if(flg)
  {
    if(!strlen(filename))
    {
      printf("You should give file name after -c option\n");
      exit(0);
    }
    if(region_array.size()>0)
    {
      printf("TIF dump to %s\n",filename);
      to_cgns(filename);
    }
    else
    {
      printf("no mesh and doping information in this TIF file.\n");
    }
  }

  PetscOptionsGetString(PETSC_NULL,"-m",filename,63,&flg);
  if(flg)
  {
    if(!strlen(filename))
    {
      printf("You should give file name after -m option\n");
      exit(0);
    }
    if(parameter_array.size()>0)
    {
      FILE *fp=fopen(filename,"w");
      int region_index=-1;
      for(size_t i=0;i<parameter_array.size();i++)
      {

        if(region_index != parameter_array[i].region-1)
        {
          region_index = parameter_array[i].region-1;
          fprintf(fp,"# %s \n",region_array[region_index].type);
        }
        for(size_t j=0;j<parameter_array[i].parameter_name_array.size();j++)
          fprintf(fp,"%s = %e\n",parameter_array[i].parameter_name_array[j].c_str(),
                  parameter_array[i].parameter_value_array[j]);
      }
      fclose(fp);
      printf("Material information save to %s\n",filename);
    }
    else
    {
      printf("no material information in this TIF file.\n");
    }
  }


  PetscOptionsHasName(PETSC_NULL,"-t",&flg);
  if(flg)
  {
    if(transol.data_name_array.size()>0)
    {
      size_t index;
      printf("There are %lu solution datas of %lu values in this TIF file\n",
             transol.data_name_array.size(),transol.data_value_array.size());
      printf("solution datas are listed below\n");
      for(size_t i=0;i<transol.data_name_array.size();i++)
      {
        printf("%lu %s  ",i,transol.data_name_array[i].c_str());
        if((i+1)%4==0) printf("\n");
      }
      do
      {
        printf("\ndump which? input index number:");
        scanf("%ld",&index);
      }
      while(index>=transol.data_name_array.size());
      printf("\nThe unit is %s \n",transol.data_unit_array[index].c_str());
      printf("Data will be written to %s.dat  \n",transol.data_name_array[index].c_str());
      sprintf(filename,"%s.dat",transol.data_name_array[index].c_str());
      FILE *fp = fopen(filename,"w");
      int  datarecord = transol.data_value_array.size()/transol.data_name_array.size();
      for(int i=0;i<datarecord;i++)
        fprintf(fp,"%e\n",transol.data_value_array[index+i*transol.data_name_array.size()]);
      fclose(fp);
      printf("%d datarecords written\n",datarecord);
    }
    else
    {
      printf("no array list data in this TIF file.\n");
    }
  }

  PetscFinalize();
  return 0;
}
