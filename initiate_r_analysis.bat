@echo off
goto :init

:header
    echo Process data from the MC Spike nextflow pipeline and then generate
    echo reports. Optionally provide the path to a ^<data_folder^> containing
    echo the output from the pipeline
    echo.
    goto :eof

:usage
    echo USAGE:
    echo   initiate_r_analysis.sh [data_folder] [-h]
    echo.
    echo Options:
    echo   -h, --help           shows this help
    echo.
    goto :eof

:init
    set "output_dir="

:parse
    if "%~1"=="" goto :main

    if /i "%~1"=="-h"         call :header & call :usage "%~2" & goto :end
    if /i "%~1"=="--help"     call :header & call :usage "%~2" & goto :end

    if not defined output_dir     set "output_dir=%~1"     & shift & goto :parse
    shift
    goto :parse

:main
    Rscript R\process_data.R %output_dir%
    Rscript generate_reports.R %output_dir%

:end
    call :cleanup
    exit /B

:cleanup
    set "output_dir="
    goto :eof
