#include <stdio.h>
#include <stdlib.h>
#define n 6

//PROTÓTIPOS DE FUNÇÕES 
float determinante(float (*mat)[n])
void mostrarResultado(float (*mat)[n])
void mostrarMatriz(float (*mat)[n])
void lerMatriz(float (*mat)[n])
void metodoGauss(float (*mat)[n], float (*L)[n], float (*U)[n])
void criaIdentidade(float (*I)[n])
void criaMatrizNula(float (*mat)[n])
void criaVetorNulo(float *vet)
void inverteMatriz(float (*L)[n], float(*U)[n], float (*resultante)[n])

/*
*
*
* Author: Lucas Amaral
*
*/ 

//ESCREVE RESULTADO NO ARQUIVO "SAIDA.txt"
void mostrarResultado(float (*mat)[n])
{
	int i, j;
	FILE *saida;
	saida = fopen("SAIDA.txt", "w");

	for(i=0; i<n; i++){
		for(j=0; j<n; j++){
			fprintf(saida, "%.2f ", mat[i][j]);
		}
		fprintf(saida, "\n");
	}
	fclose(saida);
}

//IMPRIME A MATRIZ NOS LOGS
void mostrarMatriz(float (*mat)[n])
{
    int i,j;
    for(i=0;i<n;i++){
        printf("|");
        for(j=0;j<n;j++){
        	if(mat[i][j]<0 || mat[i][j]>10)
	            printf(" %.2f", mat[i][j]);

    		else if ((mat[i][j]>0) || (mat[i][j]<10))   			
    			printf(" %.2f ", mat[i][j]);
        }
       	printf(" |\n");	
    }    
}

//LE MATRIX DO ARQUIVO "SISTEMA.txt"
void lerMatriz(float (*mat)[n])
{
	int i, j;
	FILE *entrada;
	entrada = fopen("SISTEMA.txt", "r");

	while(!feof(entrada)){
		for(i=0; i<n; i++){
			for(j=0; j<n; j++)
				fscanf(entrada, "%f", &mat[i][j]);			
		}
	}

	fclose(entrada);
}

//UTILIZA GAUSS E ENCONTRA AS MATRIZES L E U
void metodoGauss(float (*mat)[n], float (*L)[n], float (*U)[n])
{
   printf("=====================APLICANDO GAUS==================================");

    int i,j;
    //float L[n][n],U[n][n];
    for(i=0;i<n;i++){
        for(j=0;j<n;j++){
            L[i][j]=0;
            U[i][j]=mat[i][j]; 
        }
    }
    printf("\n");
    mostrarMatriz(U);

    for(j=0;j<n;j++){      
        printf("\nMultiplicando coluna %d\n",j+1);
        for(i=j+1;i<n;i++){
            if(U[i][j]==0)
            {
                printf("\n");
            }else
             {
                L[i][j]=U[i][j]/U[j][j];
                printf("M[%d][%d]=%.2f\n",i+1,j+1,L[i][j]);
                int coluna=0;
               
                for(coluna=j;coluna<n;coluna++){
               		printf("\n\n");     
                    mostrarMatriz(U);
                    U[i][coluna]=U[i][coluna]+U[j][coluna]*(-1*(L[i][j]));
               
                    printf("\n\n");    
                    mostrarMatriz(U);
                }      
            }      
            printf("\n");
        }
             
    }
    printf("\n");

    for(i=0;i<n;i++) 
    	L[i][i]=1;

   printf("=====================FIM DOS CALCULOS=============================\n");    
    
    printf("\nMatriz L\n");
    mostrarMatriz(L);
   
    printf("\nMatriz U\n");
    mostrarMatriz(U);
}


//CALCULA O DETERMINANTE DE UMA MATRIZ
float determinante(float (*mat)[n])
{
	float det=1;
	for(int i=0; i<n; i++)
		det=det*mat[i][i];
	return det;
}

//RETORNA A MATRIZ IDENTIDADE
void criaIdentidade(float (*I)[n])
{
	for(int i=0; i<n; i++){
		for(int j=0; j<n; j++){
			if(i==j)
				I[i][j]=1;
			else
				I[i][j]=0;
		}
	}
}

//CRIA UMA MATRIZ NULA
void criaMatrizNula(float (*mat)[n])
{
	for(int i=0; i<n; i++){
		for(int j=0; j<n; j++){
			mat[i][i]=0;
		}
	}
}

//CRIA VETOR NULO(USADO APENAS COMO AUXILIO)
void criaVetorNulo(float *vet)
{
	for (int i = 0; i < n; ++i)
	{
		vet[i]=0;
	}
}

/*
*
*	RESOLVE SISTEMAS LINEARES QUE RETORNAM A INVERSA
*		CALCULA L*X=I CUJO RESULTADO RESULTA NA MATRIZ Z
*		
*		CALCULA U*X=Z CUJO RESULTADO É A-¹
*
*/
void inverteMatriz(float (*L)[n], float(*U)[n], float (*resultante)[n])
{
	float identidade[n][n], auxZ[n][n], Z[n][n];
	criaIdentidade(identidade);
	criaMatrizNula(Z);
	criaMatrizNula(resultante);
	criaMatrizNula(auxZ);
	int k=0, it;


	//RESOVE O SISTEMA DE EQ LINEARES PARA Z
	criaMatrizNula(Z);
	for(int k=0; k<n; k++){
		for(int i=0; i<n; i++){
			for(int j=0, it=1; j<n; j++, it++){	
				if(i==0){
					if(i=!j)
						auxZ[k][i] = auxZ[k][i] - L[i][j];
					else if (i==j)
						auxZ[k][i] = auxZ[k][i] + identidade[i][k];
				}
				if(i>0){
					if(i!=j)
						auxZ[k][i] = auxZ[k][i] - L[i][j]*L[i][it-j];
					else if(i==j)
						auxZ[k][i] = auxZ[k][i];
				}
			}
			auxZ[k][i] = (auxZ[k][i]+ identidade[i][k])/L[i][i];
		}
	}
	printf("\n\n");
	mostrarMatriz(L);

	printf("\n\n");
	mostrarMatriz(auxZ);

	//RESOLVE O SISTEMA DE EQUAÇÕES LINEARES PARA A-¹
	for(int k=0; k<n; k++){
		for(int i=0; i<n; i++){
			for(int j=0, it=1; j<n; j++, it++){	
				if(i==0){
					if(i=!j)
						resultante[k][i] = resultante[k][i] - U[i][j];
					else if (i==j)
						resultante[k][i] = resultante[k][i] + Z[i][k];
				}
				if(i>0){
					if(i!=j)
						resultante[k][i] = resultante[k][i] - U[i][j]*U[i][it-j];
					else if(i==j)
						resultante[k][i] = resultante[k][i];
				}
			}
			resultante[k][i] = (resultante[k][i] + Z[i][k])/L[i][i];
		}
	}
	printf("\n\n");
	mostrarMatriz(resultante);


	/*printf("Vetor: \n");
	for (int i = 0; i < n; i++)
	{
		printf("%.2f\n", auxZ[i]);
	}*/

	/*
	for(int k=0; k<n; k++){
		for(int i=0; i<n; i++){
			for(int j=0; j<n; j++){
				for(int coluna=0; coluna<n; coluna++){
					if(i!=coluna)
						auxZ[i] =  identidade[i][j] + auxZ[i] - L[i][coluna]; 
						
					auxZ[i] = auxZ[i]/L[i][i];
				}
			}
		}
		for(int l=0; l<n; l++)
			Z[l][k] = auxZ[l];
	}*/
	/*
	for(int i=n-1; i>=0; i++){
		for(int j=n-1; j>=0; j++){
			for(int coluna=n-1; coluna>=0; coluna++){
				Z[i][j] = Z[i][j] + identidade[i][j] - auxZ[i][coluna];
			}
		}
	}*/
	
	/*printf("\n\nMatriz Transiente Z:\n");
	mostrarMatriz(Z);
	*/

	/*for(int i=0; i<n; i++){
		printf("%.2f\n", auxZ[i]);
	}*/
/*
	printf("\n\nMatriz Resultante Z:\n");
	mostrarMatriz(Z);*/

}

int main()
{
	//Declaração da matriz A, L e U
    float A[n][n], L[n][n], U[n][n], resultante[n][n];

    //Le matriz do arquivo SISTEMA.txt
    lerMatriz(A);

    //Mostrar A
    printf("Matriz A\n");
    mostrarMatriz(A);
    
    //Aplicar o Método de Gauss
    metodoGauss(A, L, U);

    //Tiramos o determinate a partir da matriz U
    printf("determinante= %.2f\n\n", determinante(U));
    
    //Mostrar a matriz A novamente
    printf("Matriz Original A\n");
    mostrarMatriz(A);

    inverteMatriz(L, U, resultante);

    mostrarResultado(resultante);
   
    return 0;  
}