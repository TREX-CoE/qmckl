#!/bin/bash
#
# Creates all the dependencies from the org-mode files

if [[ -z ${srcdir} ]] ; then
   echo "Error: srcdir environment variable is not defined"
   exit 1
fi


WD=$PWD

function make_src()
{

    cd ${srcdir}

    declare -A DEPS DEPS_ORG

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
    TANGLED_FILES=

    for org in org/*.org ; do
        tangled="\$(srcdir)/${org%.org}.tangled"

        i=${org%.org}
        c_test_o="\$(srcdir)/src/test_${i}.o"
        f_test_o="\$(srcdir)/src/test_${i}_f.o"
        c_test="\$(srcdir)/src/test_${i}.c"
        f_test="\$(srcdir)/src/test_${i}_f.f90"

        i="\$(srcdir)/src/${i#org/}"
        c="${i}.c"
        o="${i}.o"
        h_func="${i}_func.h"
        h_type="${i}_type.h"
        h_private_func="${i}_private_func.h"
        h_private_type="${i}_private_type.h"
        f90="${i}_f.f90"
        fo="${i}_f.o"
        fh_func="${i}_fh_func.f90"
        fh_type="${i}_fh_type.f90"

        ORG_FILES+="\$(srcdir)/$org "
        TANGLED_FILES+="$tangled "
        DEPS_ORG["\$(srcdir)/$org"]=$tangled

        grep -q "(eval c)" $org
        if [[ $? -eq 0 ]] ; then
            DEPS[$c]+=" $tangled"
            DEPS[$o]+=" $c "
            C_FILES+=" $c"
            C_O_FILES+=" $o"
        fi

        grep -q "(eval h_func)" $org
        if [[ $? -eq 0 ]] ; then
            DEPS[$o]+=" $h_func"
            DEPS[$h_func]+=" $tangled"
            H_FUNC_FILES+=" $h_func"
        fi

        grep -q "(eval h_type)" $org
        if [[ $? -eq 0 ]] ; then
            DEPS[$o]+=" $h_type"
            DEPS[$h_type]+=" $tangled"
            H_TYPE_FILES+=" $h_type"
        fi

        grep -q "(eval h_private_type)" $org
        if [[ $? -eq 0 ]] ; then
            DEPS[$o]+=" $h_private_type"
            DEPS[$h_private_type]+=" $tangled"
            H_PRIVATE_TYPE_FILES+=" $h_private_type"
        fi

        grep -q "(eval h_private_func)" $org
        if [[ $? -eq 0 ]] ; then
            DEPS[$o]+=" $h_private_func"
            DEPS[$h_private_func]+=" $tangled"
            H_PRIVATE_FUNC_FILES+=" $h_private_func"
        fi

        grep -q "(eval f)" $org
        if [[ $? -eq 0 ]] ; then
            DEPS[$f90]+="$tangled "
            DEPS[$fo]+="$f90 "
            F_FILES+=" $f90"
        fi

        grep -q "(eval fh_func)" $org
        if [[ $? -eq 0 ]] ; then
            DEPS[$fo]+=" $fh_func"
            DEPS[$fh_func]+=" $tangled"
            FH_FUNC_FILES+=" $fh_func"
        fi

        grep -q "(eval fh_type)" $org
        if [[ $? -eq 0 ]] ; then
            DEPS[$fo]+=" $fh_type"
            DEPS[$fh_type]+=" $tangled"
            FH_TYPE_FILES+=" $fh_type"
        fi

        grep -q "(eval c_test)" $org
        if [[ $? -eq 0 ]] ; then
            DEPS["$c_test"]="${tangled} "
            DEPS["$c_test_o"]+=" $c_test $o"
            C_TEST_FILES+=" $c_test"
        fi

        grep -q "(eval f_test)" $org
        if [[ $? -eq 0 ]] ; then
            DEPS["$f_test"]+="${tangled} "
            DEPS["$f_test_o"]+=" $f_test $fo"
            F_TEST_FILES+=" $f_test"
        fi
    done

    for org in org/*.org ; do
        i=${org%.org}
        i="\$(srcdir)/src/${i#org/}"
        c="${i}.c"
        o="${i}.o"
        fo="${i}_f.o"
        for i in ${!DEPS[@]} ; do
            extension="${i##*.}"
            grep -q "$i" $org
            if [[ $? -ne 0 ]] ; then
                if [[ "$extension" == h ]] ; then
                    DEPS[$o]+=" $i"
                elif [[ "$extension" == f90 ]] ; then
                    DEPS[$fo]+=" $i"
                fi
            fi
        done
    done

    for f in ${!DEPS[@]} ; do
        if [[ "${f/_f.o/_f.ox}" == ${f}x ]] ; then
            DEPS["${f}"]+=" qmckl_f.o"
        elif [[ "${f/.o/.ox}" == ${f}x ]] ; then
            DEPS["$f"]+=" \$(qmckl_h)"
        fi
    done

    OUTPUT=${WD}/generated.mk
    echo > ${OUTPUT}
    echo "## Source files" > ${OUTPUT}
    echo >> ${OUTPUT}
    echo "ORG_FILES=${ORG_FILES}" >> ${OUTPUT}
    echo "TANGLED_FILES=${TANGLED_FILES}" >> ${OUTPUT}
    echo "C_FILES=${C_FILES}" >> ${OUTPUT}
    echo "F_FILES=${F_FILES}" >> ${OUTPUT}
    echo "C_O_FILES=${C_O_FILES}" >> ${OUTPUT}
    echo "F_O_FILES=${F_O_FILES}" >> ${OUTPUT}
    echo "FH_FUNC_FILES=${FH_FUNC_FILES}" >> ${OUTPUT}
    echo "FH_TYPE_FILES=${FH_TYPE_FILES}" >> ${OUTPUT}
    echo "H_FUNC_FILES=${H_FUNC_FILES}" >> ${OUTPUT}
    echo "H_TYPE_FILES=${H_TYPE_FILES}" >> ${OUTPUT}
    echo "H_PRIVATE_FUNC_FILES=${H_PRIVATE_FUNC_FILES}" >> ${OUTPUT}
    echo "H_PRIVATE_TYPE_FILES=${H_PRIVATE_TYPE_FILES}" >> ${OUTPUT}
    echo "C_TEST_FILES=${C_TEST_FILES}" >> ${OUTPUT}
    echo "F_TEST_FILES=${F_TEST_FILES}" >> ${OUTPUT}
    echo >> ${OUTPUT}

    echo >> ${OUTPUT}
    echo "## Org-mode inherited dependencies" >> ${OUTPUT}
    echo >> ${OUTPUT}
    for f in ${!DEPS_ORG[@]} ; do
        echo ${DEPS_ORG[$f]}: $f
        echo "	\$(tangle_verbose)\$(srcdir)/tools/tangle.sh $f"
        echo ""
    done >> ${OUTPUT}
    echo >> ${OUTPUT}

    echo >> ${OUTPUT}
    echo "## Source dependencies" >> ${OUTPUT}
    echo >> ${OUTPUT}
    for f in ${!DEPS[@]} ; do
        echo ${f}: ${DEPS[$f]}
    done | sort >> ${OUTPUT}
}


make_src

