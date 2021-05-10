#!/bin/bash


function org_files() {
    echo ORG_FILES=$(echo *.org)
}

declare -A DEPS

C_FILES=
F_FILES=
FH_FUNC_FILES=
FH_TYPE_FILES=
H_FUNC_FILES=
H_TYPE_FILES=
H_PRIVATE_FUNC_FILES=
H_PRIVATE_TYPE_FILES=
C_TEST_FILES=
F_TEST_FILES=

for org in qmckl_*.org ; do
    i=${org%.org}
    c=${i}.c
    h_func=${i}_func.h
    h_type=${i}_type.h
    h_private_func=${i}_private_func.h
    h_private_type=${i}_private_type.h
    f90=${i}_f.f90
    fh_func=${i}_fh_func.f90
    fh_type=${i}_fh_type.f90

    grep -q "(eval c)" $org
    if [[ $? -eq 0 ]] ; then
        DEPS["$c"]+=" $org"
        C_FILES+=" $c"
    fi

    grep -q "(eval h_func)" $org
    if [[ $? -eq 0 ]] ; then
        DEPS[$h_func]+=" $org"
        DEPS[$c]+=" $h_func"
        H_FUNC_FILES+=" $h_func"
    fi

    grep -q "(eval h_type)" $org
    if [[ $? -eq 0 ]] ; then
        DEPS[$h_type]+=" $org"
        DEPS[$c]+=" $h_type"
        H_TYPE_FILES+=" $h_type"
    fi

    grep -q "(eval h_private_type)" $org
    if [[ $? -eq 0 ]] ; then
        DEPS[$h_private_type]+=" $org"
        DEPS[$c]+=" $h_private_type"
        H_PRIVATE_TYPE_FILES+=" $h_private_type"
    fi

    grep -q "(eval h_private_func)" $org
    if [[ $? -eq 0 ]] ; then
        DEPS[$h_private_func]+=" $org"
        DEPS[$c]+=" $h_private_func"
        H_PRIVATE_FUNC_FILES+=" $h_private_func"
    fi

    grep -q "(eval f)" $org
    if [[ $? -eq 0 ]] ; then
        DEPS[$f90]+=" $org"
        F_FILES+=" $f90"
    fi

    grep -q "(eval fh_func)" $org
    if [[ $? -eq 0 ]] ; then
        DEPS[$fh_func]+=" $org"
        DEPS[$f90]+=" $fh_func"
        FH_FUNC_FILES+=" $fh_func"
    fi

    grep -q "(eval fh_type)" $org
    if [[ $? -eq 0 ]] ; then
        DEPS[$fh_type]+=" $org"
        DEPS[$f90]+=" $fh_type"
        FH_TYPE_FILES+=" $fh_type"
    fi
done

for org in qmckl_*.org ; do
    i=${org%.org}
    c=${i}.c
    f90=${i}.f90
    c_test=test_${i}.c
    f_test=test_${i}_f.f90
    grep -q "(eval c_test)" $org
    if [[ $? -eq 0 ]] ; then
        DEPS[$c_test]+=" $org ${DEPS[$c]}"
        C_TEST_FILES+=" $c_test"
    fi

    grep -q "(eval f_test)" $org
    if [[ $? -eq 0 ]] ; then
        DEPS[$f_test]+=" $org ${DEPS[$f90]}"
        F_TEST_FILES+=" $f_test"
    fi
done

for org in ${ORG_FILES} ; do
    i=${org%.org}
    c=${i}.c
    f90=${i}.f90
    for f in ${!DEPS[@]} ; do
        extension="${f##*.}"
        grep -q "$f" $org
        if [[ $? -ne 0 ]] ; then
            if [[ extension == ".h" ]] ; then
                DEPS[$c]+=" $f"
            elif [[ extension == ".f90" ]] ; then
                DEPS[$f90]+=" $f"
            fi
        fi
    done
done

for f in ${!DEPS[@]} ; do
    if [[ "${f/_f.f90/_f.f90x}" == "${f}x" ]] ; then
        DEPS["$f"]+=" qmckl_f.o"
    fi
done

echo > generated.mk
echo "C_FILES=${C_FILES}" >> generated.mk
echo "F_FILES=${F_FILES}" >> generated.mk
echo "FH_FUNC_FILES=${FH_FUNC_FILES}" >> generated.mk
echo "FH_TYPE_FILES=${FH_TYPE_FILES}" >> generated.mk
echo "H_FUNC_FILES=${H_FUNC_FILES}" >> generated.mk
echo "H_TYPE_FILES=${H_TYPE_FILES}" >> generated.mk
echo "H_PRIVATE_FUNC_FILES=${H_PRIVATE_FUNC_FILES}" >> generated.mk
echo "H_PRIVATE_TYPE_FILES=${H_PRIVATE_TYPE_FILES}" >> generated.mk
echo "C_TEST_FILES=${C_TEST_FILES}" >> generated.mk
echo "F_TEST_FILES=${F_TEST_FILES}" >> generated.mk
echo >> generated.mk

for f in ${!DEPS[@]} ; do
    echo ${f}: ${DEPS[$f]}
done | sort >> generated.mk


