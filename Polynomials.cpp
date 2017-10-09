#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "Polynomials.h"
#include "GenericMatrix.h"
#include "GenericMatrix.cpp"
#include <QDebug>


int findmaxpowers(char * input,int * x,int * y){ //finds the max powers of x,y of the polynomial
	char check,* powerstr;
	char * tmp = input;
	int i=0;
	int counter=0;
	int power=1;
	*x=0;
	*y=0;
	int k =strlen(input);
	if((input[0]=='-')||(input[0]=='+'))
	{
  		input++;
 	}
	while(input[0]!='\n')
	{
		if(input[0]=='^'){	//is power
			counter=0;
			char * tmp=input+1;
			while((tmp[0]!='\n')&&(tmp[0]!='+')&&(tmp[0]!='*')&&(tmp[0]!='-')){	//length of the power string
				counter++;
				tmp++;
			}
			if(counter==0)
				power=1;
			else{
				powerstr=(char *)malloc(sizeof(char)*counter+1);  //+1 for \n
				strncpy(powerstr,input+1,counter);
				powerstr[counter]='\n';
				power=atoi(powerstr);
				free(powerstr);
			}
			if((input-1)[0]=='x'){	//power of x or y
				if(power>(*x))
					*x=power;
			}else if((input-1)[0]=='y'){
				if(power>(*y))
					*y=power;
			}else
				printf("error \n");
			input+=counter;		//pass processed characters
		}
		if(((input[0]=='x')||(input[0]=='y')) && (input[1]!='^')){	//x or y(not y^1 or x^1)
			if (input[0]=='x')
				if(1>(*x))
					*x=1;
			if (input[0]=='y')
				if(1>(*y))
					*y=1;
			}
		input++;
	}
}

int findnext(char * buffer,char const * n){	//finds 1->next group of the polynomial to be processed	("nextgroup")
	char * bf=buffer;						//		2->power of a variable(in string form)			("power")
	int i=0;								//      3->coefficient of a set of variables			("coeff")
	char check;								//returns the lenght of each
	if (strcmp(n,"nextgroup")==0){		
		while(bf[i]!=EOF){
			check=bf[i];
			if(check=='\n')					//end	
				return -1;
			if(check=='+' || check=='-'){	//groups are divided by -,+
				//printf("%d \n",i);
				return i;
			}
			i++;
		}
	}else if(strcmp(n,"power")==0){
		if(bf[i]==EOF)						//end of the buffer and is x or y(not x^1 or y^1)
			return -2;
		while(bf[i]!=EOF){
			check=bf[i];
			if(check=='\n')
				return -1;
			if(check=='*'){					//the power ends before a *
				//printf("%d \n",i);
				return i;
			}
			i++;
		}
	}else if(strcmp(n,"coeff")==0){
		if(bf[i]=='x'||bf[i]=='y')			// is x or y (not 1*x or 1*y)
			return 0;
		while(bf[i]!=EOF){
			check=bf[i];
			if(check=='\n')
				return -1;
			if(check=='*'){					//the coeff ends before a *
				return i;
			}
			i++;
		}
	}
}

int findpower(char * buffer,int * pass){  	//returns the number of the power(transforms it from string) , pass->number of characters the calling funtion has to pass after
	char * powerstr;
	int power;
	int powerlen=findnext(buffer,"power");	//length of the power string
	*pass=0;
	if(powerlen==0)							//y not y^1 at the end of the buffer
		return 1;
	if(powerlen==-2)						//x*y not x^1*y in the middle of the buffer
		return 1;
	if(powerlen == -1){						//power is the last to be processed in the buffer
		powerlen = strlen(buffer);
		powerstr = (char*)malloc(sizeof(char)*powerlen);
		strncpy(powerstr,buffer,powerlen);	
	}
	else{
		powerstr=(char *)malloc(sizeof(char)*powerlen+1); 	//+1 for \n
		strncpy(powerstr,buffer,powerlen);
		powerstr[powerlen]='\n';
	}
	power=atof(powerstr);
	buffer+=powerlen;
	*pass=powerlen;
	free(powerstr);
	return power;	
}

double findcoeff(char * buffer,int * pass){	//returns the number of the coeff(transforms it from string) , pass->number of characters the calling funtion has to pass after
	double coeff;
	char * coeffstr;
	int coefflen=findnext(buffer,"coeff");
	*pass=0;
	if (coefflen==0)					//x*y not 1*x*y
		return 1;
	if(coefflen == -1){						//power is the last to be processed in the buffer
		coefflen = strlen(buffer);
		coeffstr = (char*)malloc(sizeof(char)*coeff);
		strncpy(coeffstr,buffer,coefflen);	
	}
	else{
		coeffstr = (char*)malloc(sizeof(char)*coeff+1);		//+1 for \n
		strncpy(coeffstr,buffer,coefflen);	
		coeffstr[coefflen]='\n';
	}
	coeff=atof(coeffstr);
	*pass=coefflen;
	free(coeffstr);
	return coeff;	
}

void createpolmatrix(char * buffer,char * tmp_buffer,GenericMatrix<double> * pol1){
	char * hbuf;				//Used in temporary calculations
	int sizebetween,x,y,pass;	//x,y powers , pass->characters to pass , sizebetween = used to store returned values
	int sign=1;					//default sign 1(each time it is calculated again)
	if(buffer==tmp_buffer){		//the buffer to be processed starts with a + or -
		if(buffer[0]=='-'){
			sign=-1;
			buffer++;
		}
		if(buffer[0]=='+')
			buffer++;
	}
	while(buffer[0]!='\n'){		//till the end
        sizebetween=findnext(buffer,"nextgroup");	//each buffer is divided into groups consisting of coeff*variable1^power*variable2^power
        //qDebug()<<buffer;
		if(sizebetween==-1){						//last group to be processed
			sizebetween = strlen(buffer);
			hbuf = (char*)malloc(sizeof(char)*sizebetween);		//group to be processed
			char * tmp_hbuf=hbuf;
			strncpy(hbuf,buffer,sizebetween);
			double coeff=findcoeff(hbuf,&pass);					//coeff of the group		
			y=0;	
			x=0;
			if(pass != 0)										//still needs to be passed still if it is x*y (not 1*x*y)
			hbuf+=pass+1;
			if(hbuf[0]=='x'){									//first there is the x variable
			hbuf++;
			if(hbuf[0]=='^'){									//find its power
				hbuf++;
				x=findpower(hbuf,&pass);
				hbuf+=pass+1;
			}else{
				x=1;
				hbuf++;
			}
			if(hbuf[0]=='y'){									//then there is the y variable
				hbuf++;
				if(hbuf[0]=='^'){
					hbuf++;
					y=findpower(hbuf,&pass);
				}else{
					y=1;
					hbuf++;
				}
			}
		}else if(hbuf[0]=='y'){									//first there is the y variable
			hbuf++;
			if(hbuf[0]=='^'){
				hbuf++;
				y=findpower(hbuf,&pass);
				hbuf+=pass+1;
			}else{
				y=1;
				hbuf++;
			}
			if(hbuf[0]=='x'){									//then there is the y variable
				hbuf++;
				if(hbuf[0]=='^'){
					hbuf++;
					x=findpower(hbuf,&pass);
					hbuf+=pass+1;
				}else{
					x=1;
					hbuf++;
				}
			}
		}
		if ((buffer-1)[0]=='-')						//find the sign of the group
			sign=-1;
		else 
			sign=1;
		buffer=buffer+(sizebetween+1);					
		pol1->SetMatrix(y,x,sign*coeff);						//insert to the matrix
		free(tmp_hbuf);	
			break;
		}
		hbuf = (char*)malloc(sizeof(char)*sizebetween+1);		//not last group to be processed
		char * tmp_hbuf=hbuf;									//same part as before
		strncpy(hbuf,buffer,sizebetween);
		hbuf[sizebetween]='\n';
		double coeff=findcoeff(hbuf,&pass);
		y=0;
		x=0;
		if(pass != 0)
			hbuf+=pass+1;
		if(hbuf[0]=='x'){
			hbuf++;
			if(hbuf[0]=='^'){
				hbuf++;
				x=findpower(hbuf,&pass);
				hbuf+=pass+1;
			}else{
				x=1;
				hbuf++;
			}
			if(hbuf[0]=='y'){
				hbuf++;
				if(hbuf[0]=='^'){
					hbuf++;
					y=findpower(hbuf,&pass);
				}else{
					y=1;
					hbuf++;
				}
			}
		}else if(hbuf[0]=='y'){
			hbuf++;
			if(hbuf[0]=='^'){
				hbuf++;
				y=findpower(hbuf,&pass);
				hbuf+=pass+1;
			}else{
				y=1;
				hbuf++;
			}
			if(hbuf[0]=='x'){
				hbuf++;
				if(hbuf[0]=='^'){
					hbuf++;
					x=findpower(hbuf,&pass);
					hbuf+=pass+1;
				}else{
					x=1;
					hbuf++;
				}
			}
		}
		pol1->SetMatrix(y,x,sign*coeff);
		if ((buffer+sizebetween)[0]=='-')
			sign=-1;
		else 
			sign=1;
		buffer=buffer+(sizebetween+1);			//set for next group
		free(tmp_hbuf);	
	}
}

void RandomizePolynomial(GenericMatrix<double> * pol,int r,int c)
{
	int i,j,sign;
	for(i = 0 ; i<=r; i++)
	{
		for(j = 0 ; j<=c; j++)
		{
			sign = rand()%10;
			if(sign>4)
				sign=-1;
			else
				sign=1;
			pol->SetMatrix(i,j,sign*rand()%10);
		}
	}
}
