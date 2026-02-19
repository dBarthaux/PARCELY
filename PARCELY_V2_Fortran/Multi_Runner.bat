@echo off
setlocal enabledelayedexpansion

rem Set the batch file's directory
set "batch_file_directory=%~dp0"

rem Prompt user for number of runs
set /p "DataCols=Enter the number of data columns in the InputTable.csv file: "
set "DataCols=!DataCols!"

rem For each data column
for /l %%i in (1, 1, %DataCols%-1) do (
	rem Set the argument to be passed to the Python script
	set "argument1=%%i"
	rem Run the input file changer
	python Input_Changer.py !argument1!
	rem Change directory where the application file is
	cd "PARCELY_V2"
	rem Run PARCELY
	echo.|PARCELY_V2
	rem Create a new directory for the ouput files
	cd "!batch_file_directory!/PARCELY_Output"
	set "foldName=Output_%%i"
	mkdir "!foldName!"
	rem Move all the files into the new output folder
	move "*" "!foldName!"
	cd "!batch_file_directory!"
	rem pause
)
endlocal