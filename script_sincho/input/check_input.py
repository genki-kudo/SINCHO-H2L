from Bio.PDB import *
import argparse
import os
import logging
from script_sincho.basic.basic_func import *


def setting():

    input_list=[0 for i in range(11)]

    parser = argparse.ArgumentParser(prog='sincho')
    #input file setting
    parser.add_argument('-c', '--clusterdir', nargs=1, default=['./p2c_output/cluster/'], help='specify p2c output directory with pocket clusters if neccesary (default: ./p2c_output/cluster/)') 
    parser.add_argument('-o', '--outputdir', nargs=1, default=['./sincho-output/'], help='specify output directory if neccesary (default: ./sincho-output)')
    parser.add_argument('-f', '--fpocketdir', nargs=1, default=['./asphere_output/'], help='specify fpocket output directory if neccesary (default: ./asphere_output)')
    parser.add_argument('-l', '--ligand', nargs=1, required=True, help='specify ligand pdb file. ')
    parser.add_argument('-p', '--protein', nargs=1, required=True, help='specify protein pdb file. ')
    #####to me ... make -noh option about (B) and (C) #####
    #parser.add_argument('-noh', '--nohydrogen', help='you use this flag if ligand.pdb file does not include hydrogen', action='store_true')
    #####to me ... make -nop2c option about (A) #####
    #parser.add_argument('-skipp2c', '--skipp2c', help='you use this flag if you does not execute Pocket-to-Concavity', action='store_true')
    ## if you already executed P2C before SINCHO, you have this flag. 
    parser.add_argument('-log', '--logfilename', nargs=1, default=['./sincho.log'], help='specify logfile name (default: sincho.log)')
    #parser.add_argument('-w', '--weight', nargs=1, default=['0.5'], help='specify 0.0 to 1.0. default: 0.5 (no weighted in the extend score)')
    parser.add_argument('-n', '--number_cand', nargs=1, default=['10'], help='specify number of candidates. default: 10')
    #parser.add_argument('-mwsd','--mwerror',nargs=1, default=['-50'], help='correct term of mw-estimation. default: -50')

    args = parser.parse_args()

    input_list[0] = str(args.clusterdir[0])#pocket cluster dir
    input_list[1] = str(args.outputdir[0])#output dir
    input_list[2] = str(args.fpocketdir[0])#fpocket output dir
    input_list[3] = str(args.ligand[0])#ligand pdb
    input_list[4] = str(args.protein[0])#protein pdb
    #input_list[5] = str(args.nohydrogen)#no hydrogen flag
    #input_list[6] = str(args.skipp2c)#no p2c flag
    input_list[7] = str(args.logfilename[0])#log file name
    #input_list[8] = str(args.weight[0])#weight in the extend score
    input_list[9] = str(args.number_cand[0])# number of pocket selected at (B)
    #input_list[10] = str(args.mwerror[0])

    logger = logging.getLogger(__name__)
    logger.setLevel(logging.DEBUG)
    formatter = logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
    t_file(input_list[7])
    file_handler = logging.FileHandler(input_list[7])
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)

    return input_list, logger

def checks(input_list,logger):
    is_file = os.path.isfile(input_list[3])
    if is_file:
        logger.info('Ligand file>> '+input_list[3]+' >>OK.')
    else:
        return
    #if input_list[5]=='True':
        #logger.warning('This ligand file ('+input_list[3]+') does not include hydrogen atoms.(flag: -noh)')
    #if input_list[6]=='True':
        #logger.info('if you already executed P2C before SINCHO, you have this flag.')
    #else:
    logger.info('Pocket cluster directory>> '+input_list[0])
    logger.info('Fpocket output directory>> '+input_list[2])
    logger.info('Output directory>> '+input_list[1])
    if os.path.exists(input_list[1]):
        logger.info('Output directory>> '+input_list[1])
        os.makedirs(input_list[1]+"/pocket_environment/", exist_ok=True)
    else:
        logger.warning('Output directory: '+input_list[1]+' will be generated.')
        os.makedirs(input_list[1], exist_ok=True)
        os.makedirs(input_list[1]+"/pocket_environment/", exist_ok=True)

    return







