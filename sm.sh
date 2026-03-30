#!/bin/bash
#
#	Autores:	 Jorge De La Peña García 	(  jorge.dlpg@gmail.com  )
#				 Horacio Pérez Sánchez   	(  hperez@ucam.edu       )
#				 Ricardo Rodriguez schmidt	(  rrschmidt@ucam.edu    )	
#
#_______________________________________________________________________________________________________________________________

error=0		  
paraMError="" 	  
txtErrror=""  	  
informe=""	  
fecha=$(date +"%Y-%m-%d")
pathSL="GROMACS/login_node/"
path_gromacs="${PWD}/GROMACS/"

function f_help(){
	source ${pathSL}help.sh $1
	if [ ! -z "$txtErrror" ];then 
		echo ""
		echo -e "${RED}ERROR: ${BROWN}"$txtErrror"${NONE}"
		echo ""
	fi
	exit
}


#______________________________________________________________________
#
#	Funcion para indicar si existe algun error de los datos de entrada
#____________________________________________________________________
function f_Error()
{
	
	laError="GROMACS.sh,"
	case $error in
		1)  txtErrror="$laError Programa de docking incorrecto"   ; f_help;;
		2)  txtErrror="$laError Tecnica de docking incorrecto para ese programa"; f_help;;
		3)  txtErrror="$laError Extension no valida para el programa de docking"; f_help;;
		4)  txtErrror="$laError Fichero ligando o de proteina vacio"; f_help;;
		5)  txtErrror="$laError Hostanme del cluster no encontrado"; f_help;;
		6)  txtErrror="$laError No encontrada opciond de docking"; f_help;;
		7)  txtErrror="$laError numJobs, x, y, o z estan vacios"; f_help;;
		8)  txtErrror="$laError Fichero de la proteina, carpeta de ligandos o fichero de ligando NO EXISTE";f_help;;
		9)  txtErrror="$laError No se pueden usar guiones medios (-) en los nombres de ficheros de proteinas o ligandos";f_help;;
		10) txtErrror="$laError Gestor de colas erroneo, ";f_help;;
		11) txtErrror="$laError Parametro renice errorne, revise la aydua";f_help;;
		12) txtErrror="$laError Parametro Desconocido $paraMError, revise la aydua";f_help;;
	esac
}
#___________________________________________________________________________________________________________________________________________________________#
#																																							#
#															MAIN																							#
#___________________________________________________________________________________________________________________________________________________________#
#
#	parametros, colores y configuacion global
#
source ${pathSL}colors.sh			  			#colores del script
source ${pathSL}parameters.sh 		  				#parametros 
source ${pathSL}debug.sh 						#Para que funcione el debug
source ${pathSL}create_conf.sh						#crea oi carga configuracion global

source ${path_login_node}/folders_experiment.sh 	#variables folders de la prueba
#
#	Si se ejcuta con -k matara los procesos del usuario que esten activos con GROMACS
#	Util cuando se ejcuta en secuencial y se sale con control z
#
if [ "$kill_sm" != "N/A" ];then
	ps fux | grep -F "GROMACS/login_node" | awk '{print $2}' | xargs kill -9 2> /dev/null
	exit
elif [ "$mode_test" != "N/A" ];then
	if [ "$option" != "" ];then
		if [ -f "${pathSL}/test_params/test_params_${option}.sh" ]; then
			source ${pathSL}/test_params/test_params_${option}.sh
		else
			txtErrror="No parameters found for the test ${option}"
			f_help
			
		fi
	else
		source ${pathSL}/test_params/test_params_AD.sh
	fi
	
fi

source ${path_login_node}check_renice.sh 		#busca el cluster donde se encuentra (Se podria mejorar)
#echo $error

source ${path_login_node}verify_input_data.sh 	#comprueba los datos de entrada con los ficheros SW.txt y OP.txt
#echo $error

f_Error 			#se comprueba a ver si ha dado algun error
source ${path_login_node}create_dirs.sh 		#crea un directorio
#echo $error
source ${path_login_node}create_params.sh  #crea un fichero con los parametros especiales para ese programa
f_Error 								#se comprueba a ver si ha dado algun error

source ${path_login_node}upload_excel/update_excel.sh #Sube al eesl de google si se requiere#Sube al eesl de google si se requiere

if [ -f "${path_login_node}techniques/SLTechnique${option}.sh" ];then 
	source ${path_login_node}techniques/SLTechnique${option}.sh #ejecucion de la tecnica
else
	error=2
	f_Error 								#se comprueba a ver si ha dado algun error
fi











