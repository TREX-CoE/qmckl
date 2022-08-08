#!/usr/bin/env python3

import os
import sys

qmckl_h = sys.argv[1]

collect = False
process = False
get_name = False
block = []
res_str = ''
func_name = ''
arrays = {}
numbers = {}
qmckl_public_api = []
qmckl_errors = []

with open(qmckl_h, 'r') as f_in:
    for line in f_in:

        # get the errors but without the type cast because SWIG does not recognize it
        if '#define' in line and 'qmckl_exit_code' in line:
            qmckl_errors.append(line.strip().replace('(qmckl_exit_code)',''))
            continue

        if get_name:
            words = line.strip().split()
            if '(' in words[0]:
                func_name = words[0].split('(')[0]
            else:
                func_name = words[0]
            if 'get' in func_name or 'set' in func_name:
                qmckl_public_api.append(func_name)

            get_name = False

        if 'qmckl_exit_code' in line:
            words = line.strip().split()
            if len(words) > 1 and 'qmckl_exit_code' in words[0]:
                # this means that the function name is on the same line as `qmckl_exit_code`
                func_name = words[1].split('(')[0]
                if 'get' in func_name or 'set' in func_name:
                    qmckl_public_api.append(func_name)
            elif len(words) == 1:
                # this means that the function name is the first element on the next line
                get_name = True
                #continue # do not `continue` here otherwise collect is not True for some functions

            # process functions - oneliners (for arrays)
            if 'size_max' in line and ';' in line:

                tmp_list = line.split(',')
                for i,s in enumerate(tmp_list):
                    if 'size_max' in s:
                        end_str  = tmp_list[i].replace(';','').replace('\n','')
                        pattern  = f"({tmp_list[i-1]} ,{end_str}"
                        datatype = tmp_list[i-1].replace('const','').replace('*','').split()[0]
                        arrays[func_name] = {
                                'datatype' : datatype,
                                'pattern'  : pattern
                                }
                        #if 'qmckl_get_jastrow_type_nucl_vector' in func_name:
                        #    print(line)
                        #    print(pattern)
                continue

            # if size_max is not provided then the function should deal with numbers or string
            #elif 'num' in line and 'get' in func_name:
            elif ';' in line and 'get' in func_name:
                # special case
                if 'size_max' in line:
                    continue

                #print(line)

                tmp_str = line.split(',')[-1].strip()

                pattern = tmp_str.replace(')','').replace(';','')
                datatype = pattern.replace('const','').replace('*','').split()[0]

                numbers[func_name] = {
                    'datatype' : datatype,
                    'pattern'  : pattern
                    }
                continue
            # for multilne functions - append line by line to the list
            else:
                block.append(line)
                collect = True
                continue

        # if size_max is encountered within the multiline function
        if 'size_max' in line and collect:
            #if 'qmckl_get_electron_rescale_factor_en' in func_name:
            #    print("LOL")

            # this will not work for 2-line functions where array argument is on the same line as
            # func name and size_max argument is on the next line
            if not 'qmckl_exit_code' in block[-1] and not '*/' in line:
                pattern = '(' + block[-1].strip() + line.strip().replace(';','')
                datatype = pattern.replace('const','').replace('*','').replace('(','').split()[0]

            collect = False
            block = []
            arrays[func_name] = {
                    'datatype' : datatype,
                    'pattern'  : pattern
                    }
            continue

        #if 'num' in line and 'get' in func_name and not 'qmckl_get' in line and collect:
        if 'get' in func_name and not 'qmckl_get' in line and collect and ';' in line:
            #print(func_name)
            #print(line)
            pattern = line.replace(';','').replace(')','').strip()
            datatype = pattern.replace('const','').replace('*','').split()[0]

            collect = False
            block = []
            numbers[func_name] = {
                    'datatype' : datatype,
                    'pattern'  : pattern
                    }
            continue

        # stop/continue multiline function analyzer
        if collect and ')' in line:
            collect = False
            block = []
            continue
        else:
            block.append(line)
            continue


# remove buggy qmckl_get_electron_rescale_factor_en key
#arrays.pop('qmckl_get_electron_rescale_factor_en')

processed = list(arrays.keys()) + list(numbers.keys())

#for pub_func in qmckl_public_api:
    #if pub_func not in processed and 'set' not in pub_func:
        #print("TODO", pub_func)
    #print(v['datatype'])

#for k,v in numbers.items():
#    print(v)


with open("qmckl_include.i", 'w') as f_out:

    # write the list of errors as constants without the type cast
    for e in qmckl_errors:
        line = e.replace('#define', '%constant qmckl_exit_code').replace('(','=').replace(')',';')
        f_out.write(line + '\n')

    swig_type = ''
    for v in numbers.values():

        if 'int' in v['datatype']:
            swig_type = 'int'
        elif 'float' in v['datatype'] or 'double' in v['datatype']:
            swig_type = 'float'
        elif 'char' in v['datatype'] or 'bool' in v['datatype']:
            #print('SWIG, skipping', v['datatype'], v['pattern'])
            continue
        else:
            raise TypeError(f"Unknown datatype for swig conversion: {v['datatype']}")

        f_out.write(f"%apply {swig_type} *OUTPUT {{ {v['pattern']} }};\n")

    for k,v in arrays.items():
        if 'char' in v['datatype']:
            #print("String type", k, v)
            pass

        if len(v['pattern'].split(',')) != 2:
            print('Problemo', k, v)
            continue

        if 'get' in k:
            f_out.write(f"%apply ( {v['datatype']}* ARGOUT_ARRAY1 , int64_t DIM1 ) {{ {v['pattern']} }};\n")
        elif 'set' in k:
            f_out.write(f"%apply ( {v['datatype']}* IN_ARRAY1 , int64_t DIM1 ) {{ {v['pattern']} }};\n")
        #else:
            #print("HOW-TO ?", k)

