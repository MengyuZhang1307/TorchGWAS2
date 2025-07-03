#pragma once
#include "ReadBGEN.h"


/**************************************
This function is revised based on the Parse function in BOLT-LMM v2.3 source code
*************************************/

void Bgen::processBgenHeaderBlock(std::string bgenfile) 
{
    char genofile[300];
    strcpy(genofile, bgenfile.c_str());
    std::string genopath(genofile);
    if (genopath.substr(genopath.length() - 5, 5) != ".bgen") 
    {
        std::cout << "\nERROR: " << genopath << " does not have a .bgen extension. \n\n";
        exit(1);
    }

    fin = fopen(genofile, "rb");
    if (fin == 0) 
    {
        std::cerr << "\nERROR: BGEN file could not be opened.\n\n";
        exit(1);
    }

    std::cout << "General information of BGEN file: \n";
    if (!fread(&offset, 4, 1, fin)) 
    {
        std::cerr << "\nERROR: Cannot read BGEN header (offset).\n\n";
        exit(1);
    }

    uint L_H; 
    if (!fread(&L_H, 4, 1, fin)) 
    {
        std::cerr << "\nERROR: Cannot read BGEN header (LH).\n\n";
        exit(1);
    }

    if (fread(&Mbgen, 4, 1, fin)) 
    { 
        if (Mbgen <= 0) {
            std::cerr << "\nERROR: The number of variants in the BGEN file is 0.\n\n";
            exit(1);
        }
        std::cout << "Number of variants: " << Mbgen << '\n';
    } 
    else 
    { 
        std::cerr << "\nERROR: Cannot read BGEN header (M).\n\n"; 
        exit(1);
    }

    if (fread(&Nbgen, 4, 1, fin)) 
    {
        if (Nbgen <= 0) 
        {
            std::cerr << "\nERROR: The number of samples in the BGEN file is 0.\n\n";
            exit(1);
        }
        std::cout << "Number of samples: " << Nbgen << '\n';
    }
    else 
    {
        std::cerr << "\nERROR: Cannot read BGEN header (N). \n\n";
        exit(1);
    }


    char magic[5]; 
    if (!fread(magic, 1, 4, fin)) 
    { 
        std::cerr << "\nERROR: Cannot read BGEN header (magic bytes). \n\n"; 
        exit(1);
    }
    magic[4] = '\0';
    if (!(magic[0] == 'b' && magic[1] == 'g' && magic[2] == 'e' && magic[3] == 'n')) 
    {
        std::cerr << "\nERROR: BGEN file's four magic number bytes does not match 'b' 'g' 'e' 'n'.\n\n";
        exit(1);
    }

    fseek(fin, L_H - 20, SEEK_CUR);         

    uint flags; 
    if (!fread(&flags, 4, 1, fin)) 
    {
        std::cerr << "\nERROR: Cannot read BGEN header (flags). \n\n";
        exit(1); 
    }


    // The header block - flag definitions
    CompressedSNPBlocks = flags & 3;
    
    switch (CompressedSNPBlocks)
    {
    case 0:
        std::cout << "Genotype Block Compression Type: Uncompressed\n";
        break;
    case 1:
        std::cout << "Genotype Block Compression Type: Zlib\n";
        break;
    case 2:
        std::cout << "Genotype Block Compression Type: Zstd\n";
        break;
    default:
        std::cout << "\nERROR: BGEN compression flag must be 0 (uncompressed), 1 (zlib compression), or 2 (zstd compression). Value in file: " << CompressedSNPBlocks << ".\n\n";
        exit(1);
    }

    Layout = (flags >> 2) & 0xf; 
    std::cout << "Layout: " << Layout << '\n';
    if (Layout != 1U && Layout != 2U) 
    {
        std::cerr << "\nERROR: BGEN layout flag must be 1 or 2.\n\n";
        exit(1);
    }

    SampleIdentifiers = flags >> 31; 
    std::cout << "Sample Identifiers Present: ";  
    SampleIdentifiers == 0 ? std::cout << "False \n" : std::cout << "True \n";
    if (SampleIdentifiers != 0 && SampleIdentifiers != 1) 
    {
        std::cerr << "\nERROR: BGEN sample identifier flag must be 0 or 1.\n\n";
        exit(1);
    }
    
}



/**********************************************************************************
This function is revised based on the Parse function in BOLT-LMM v2.3 source code
***********************************************************************************/

// This functions reads the sample block of BGEN v1.1, v1.2, and v1.3. Also finds which samples to remove if they have missing values in the pheno file.
void Bgen::processBgenSampleBlock(Bgen bgen, char samplefile[300], bool useSample, UMap_str_VV_string phenomap, std::string phenoMissingKey, int numSelCol, int samSize) 
{
    int k = 0;
    std::unordered_set<int> genoUnMatchID;
    std::vector<std::string> tempID;
    new_phenodata.resize(samSize);
    std::vector<double> new_covdata_orig(samSize * (numSelCol+1));
    if ((bgen.SampleIdentifiers == 0) || useSample) {
        if (bgen.SampleIdentifiers == 0 && !useSample) {
            std::cerr << "\nERROR: BGEN file does not contain sample identifiers. A .sample file is required. \n"
                << "       See https://www.well.ox.ac.uk/~gav/qctool/documentation/sample_file_formats.html for .sample file format. \n\n";
            exit(1);
        }

        std::ifstream fIDMat;
        fIDMat.open(samplefile);
        if (!fIDMat.is_open()) {
            std::cerr << "\nERROR: Sample file could not be opened.\n\n";
            exit(1);
        }

        std::string IDline;
        std::getline(fIDMat, IDline);
        std::getline(fIDMat, IDline);
        uint nSamples = 0;
        while (getline(fIDMat, IDline)) {
            nSamples++;
        }
        if (nSamples != bgen.Nbgen) {
            std::cout << "\nERROR: Number of sample identifiers in .sample file (" << nSamples << ") does not match the number of samples specified in BGEN file (" << bgen.Nbgen << ").\n\n";
            exit(1);
        }
        else {
            fIDMat.clear();
            fIDMat.seekg(0, fIDMat.beg);
            getline(fIDMat, IDline);
            getline(fIDMat, IDline);
        }

        for (uint m = 0; m < bgen.Nbgen; m++) {
            // IDMatching
            getline(fIDMat, IDline);
            std::istringstream iss(IDline);
            std::string strtmp;
            iss >> strtmp;
            //AllsampleIDs before matching
            sampleID_all.push_back(strtmp);

            int itmp = k;
            
            if (phenomap.find(strtmp) != phenomap.end()) 
            {
                auto& tmp_valvecs = phenomap[strtmp]; 

                for (const auto& tmp_valvec : tmp_valvecs) {
                    // Check for missing phenotype values in the current vector
                    if (find(tmp_valvec.begin(), tmp_valvec.end(), phenoMissingKey) == tmp_valvec.end() &&
                        find(tmp_valvec.begin(), tmp_valvec.end(), "") == tmp_valvec.end()) 
                    {     
                        sscanf(tmp_valvec[0].c_str(), "%lf", &new_phenodata[k]);
                        new_covdata_orig[k * (numSelCol + 1)] = 1.0;
                        for (int c = 0; c < numSelCol; c++) 
                        {
                            sscanf(tmp_valvec[c + 1].c_str(), "%lf", &new_covdata_orig[k * (numSelCol + 1) + c + 1]);
                        }
                        sampleID.push_back(strtmp);
                        k++;
                    }
                }
            }
            // save the index with unmatched ID into genoUnMatchID.
            if (itmp == k) 
            {
                genoUnMatchID.insert(m);
            }
        } 
        fIDMat.close();
    }

    if ((bgen.SampleIdentifiers == 1) && !useSample) {

        uint maxLA = 65536;
        char* samID = new char[maxLA + 1];

        uint LS1;  
        if (!fread(&LS1, 4, 1, bgen.fin)) {
            std::cerr << "\nERROR: Cannot read BGEN sample block (LS).\n\n";
            exit(1);
        }

        uint Nrow; 
        if (!fread(&Nrow, 4, 1, bgen.fin)) {
            std::cerr << "\nERROR: Cannot read BGEN sample block (N).\n\n";
            exit(1);
        }
        if (Nrow != bgen.Nbgen) {
            std::cerr << "\nERROR: Number of sample identifiers (" << Nrow << ") does not match number of samples specified in BGEN file (" << bgen.Nbgen << ").\n\n";
            exit(1);
        }


        for (uint m = 0; m < bgen.Nbgen; m++) {
            ushort LSID; 
            if (!fread(&LSID, 2, 1, bgen.fin)) { 
                std::cerr << "\nERROR: Cannot read BGEN sample block (LSID).\n\n"; 
                exit(1); 
            }
            if (!fread(samID, 1, LSID, bgen.fin)) {
                std::cerr << "\nERROR: Cannot read BGEN sample block (sample id).\n\n";
                exit(1);
            }
            samID[LSID] = '\0';

            std::string strtmp(samID);
            sampleID_all.push_back(strtmp);
            int itmp = k;
            
            if (m < 5) {
                tempID.push_back(strtmp);
            }

            if (phenomap.find(strtmp) != phenomap.end()) 
            {
                auto& tmp_valvecs = phenomap[strtmp]; 

                for (const auto& tmp_valvec : tmp_valvecs) {
                    // Check for missing phenotype values in the current vector
                    if (find(tmp_valvec.begin(), tmp_valvec.end(), phenoMissingKey) == tmp_valvec.end() &&
                        find(tmp_valvec.begin(), tmp_valvec.end(), "") == tmp_valvec.end()) 
                    {     
                        sscanf(tmp_valvec[0].c_str(), "%lf", &new_phenodata[k]);
                        new_covdata_orig[k * (numSelCol + 1)] = 1.0;
                        for (int c = 0; c < numSelCol; c++) 
                        {
                            sscanf(tmp_valvec[c + 1].c_str(), "%lf", &new_covdata_orig[k * (numSelCol + 1) + c + 1]);
                        }
                        sampleID.push_back(strtmp);
                        k++;
                    }
                }
            }

            if (itmp == k) 
            {
                genoUnMatchID.insert(m);
            }
        }
        delete[] samID;
    } // end SampleIdentifiers == 1


    // After IDMatching, resizing phenodata and covdata, and updating samSize;

    new_phenodata.resize(k);
    new_covdata_orig.resize(k * (numSelCol + 1));
    samSize = k;

    if (samSize == 0) 
    {
        std::cerr << "\nERROR: Sample size changed from " << samSize + genoUnMatchID.size() << " to " << samSize << ".\n\n";
        if (bgen.SampleIdentifiers == 1 && !useSample) {
            int print_i = 5;
            if (bgen.Nbgen < 5)
            { 
                print_i = bgen.Nbgen; 
            }
            std::cout << "ID matching was done using the BGEN sample identifier block. \nHere are the first " << print_i << " sample identifiers in BGEN file: \n";
            for (int i = 0; i < print_i; i++) 
            {
                std::cout << " " << tempID[i] << "\n";
            }
        }

        if (bgen.SampleIdentifiers == 0 || useSample) 
        {
            std::cout << "Check if sample IDs are consistent between the phenotype file and sample file, or check if (--sampleid-name) is specified correctly. \n\n";
        }
        exit(1);
    }


    int ii = 0;
    include_idx.resize(samSize);
    for (uint i = 0; i < bgen.Nbgen; i++) {
        if (genoUnMatchID.find(i) == genoUnMatchID.end()) 
        {
             include_idx[ii] = i;
              ii++;
        }
    }

    std::cout << "****************************************************************************\n";
    if (genoUnMatchID.empty()) 
    {
        std::cout << "After processes of sample IDMatching and checking missing values, the sample size does not change.\n\n";
    }
    else 
    {
        std::cout << "After processes of sample IDMatching and checking missing values, the sample size changes from "
            << samSize + genoUnMatchID.size() << " to " << samSize << ".\n\n";
    }
    std::cout << "Sample IDMatching and checking missing values processes have been completed.\n";
    std:: cout << "New pheno and covariate data vectors with the same order of sample ID sequence of geno data are updated.\n";
    std::cout << "****************************************************************************\n";


    new_samSize = samSize;
    if (new_samSize<(numSelCol+1) || new_samSize == (numSelCol+1))
    {
        std::cout << "\nERROR: The sample size should be greater than the number of predictors!" <<std::endl;
        exit(1);
    }

    // //the first column of matcovX is Y
    // MatrixXd matcovX (samSize,(numSelCol+1));
    // for (int i=0; i<samSize; i++){    
    //     for (int j=0; j<(numSelCol+1); j++) {
    //       matcovX(i,j) =new_covdata_orig [i * (numSelCol+1) +j];
    //     }
    // }
    // Eigen::HouseholderQR<MatrixXd> qr;
    // qr.compute(matcovX);
    // Eigen::MatrixXd R = qr.matrixQR();
    // int colR=R.cols();
    // VectorXd diagR (colR);
    // for (int i=0; i<colR; i++){
    //     diagR(i)=abs(R(i,i));
    // }

    // double sqrtEps =sqrt(std::numeric_limits<double>::epsilon());
    // double maxdiag = *std::max_element( diagR.begin(), diagR.end() ) ;
    // double colinear_cut = abs(maxdiag * sqrtEps);
    // for (int i=0; i<colR; i++){
    //     if (abs(diagR(i)) < colinear_cut){
    //         excludeCol.push_back(i);    
    //     }
    // }
    // matcovX.resize(0,0);
    // R.resize(0,0);

    // int NumExcludeCol = excludeCol.size();
    // if (excludeCol.size()>0){        
    //     vector <int> remove_colinear;
    //     for (int i=0; i<excludeCol.size(); i++){
    //         for (int j=0; j<samSize; j++) {
    //             remove_colinear.push_back(j * (numSelCol+1) + excludeCol[i]);
    //         }
    //     }

    //     numSelCol=numSelCol- excludeCol.size();
    //     new_covdata.resize(samSize * (numSelCol+1));
    //     vector<double> temp;
    //     for (int i=0; i<new_covdata_orig.size(); i++)
    //     {
    //         if (std::find(remove_colinear.begin(), remove_colinear.end(), i) == remove_colinear.end())
    //         {
    //             temp.push_back(new_covdata_orig[i]);
                
    //         }
    //     }
    //     new_covdata = temp;
    // } 
    // else 
    // {
    //         new_covdata.resize(samSize * (numSelCol+1));
    //         new_covdata = new_covdata_orig;
    // }

}



/***********************************************************************************
This function contains code that is revised based on BOLT-LMM v2.3 source code
************************************************************************************/

// This function reads just the variant block for BGEN files version v1.1, v1.2, and v1.3 and is used to grab the byte where the variant begins.
//    Necesary when there's no bgen index file.
void Bgen::getPositionOfBgenVariant(Bgen bgen, int threads, std::string includeVariantFile, bool doFilters) {


    int count = 0;
    uint CompressedSNPBlocks = bgen.CompressedSNPBlocks;
    uint Layout = bgen.Layout;
    uint offset = bgen.offset;
    uint Mbgen = bgen.Mbgen;
    uint nSNPS = Mbgen;
    uint Nbgen = bgen.Nbgen;
    uint maxLA = 65536;
    char* snpID   = new char[maxLA + 1];
    char* rsID    = new char[maxLA + 1];
    char* chrStr  = new char[maxLA + 1];
    char* allele1 = new char[maxLA + 1];
    char* allele0 = new char[maxLA + 1];
    std::string IDline;

    extTypes::StringSet includeVariant;
    std::vector<std::vector<uint>> includeVariantIndex;
    bool checkSNPID = false;
    bool checkRSID = false;
    bool checkInclude = false;
    int ret;

    if (doFilters) 
    {

        filterVariants = true;
        if (!includeVariantFile.empty()) 
        {
            checkInclude = true;

            std::ifstream fInclude;
            fInclude.open(includeVariantFile);
            if (!fInclude.is_open()) {
                std::cerr << "\nERROR: The file (" << includeVariantFile << ") could not be opened.\n\n";
                exit(1);
            }

            std::string vars;
            getline(fInclude, IDline);
            std::transform(IDline.begin(), IDline.end(), IDline.begin(), ::tolower);
            IDline.erase(std::remove(IDline.begin(), IDline.end(), '\r'), IDline.end());
            if (IDline == "snpid") {
                std::cout << "An include snp file was detected... \nIncluding SNPs for analysis based on their snpid... \n";
                checkSNPID = true;
            }
            else if (IDline == "rsid") {
                std::cout << "An include snp file was detected... \nIncluding SNPs for analysis based on their rsid... \n";
                checkRSID = true;
            }
            else {
                std::cerr << "\nERROR: Header name of " << includeVariantFile << " must be snpid or rsid.\n\n";
                exit(1);
            }

            while (fInclude >> vars) 
            {
                if (includeVariant.find(vars) != includeVariant.end()) 
                {
                    std::cout << "\nERROR: " << vars << " is a duplicate variant in " << includeVariantFile << ".\n\n";
                    exit(1);
                }
                includeVariant.insert(vars);
                count++;
            }
            nSNPS = count;
            std::cout << "Detected " << nSNPS << " variants to be used for analysis... \nAll other variants will be excluded.\n\n\n";
            std::cout << "Detected" << std::thread::hardware_concurrency() << " available thread(s)...\n"; 
            if (nSNPS < threads) {
                std::cout << "Number of variants (" << nSNPS << ") is less than the number of specified threads (" << threads << ")...\n";
                threads = nSNPS;
                std::cout << "Using " << threads << " for multithreading... \n\n";
            }
            else {
                std::cout << "Using " << threads << " for multithreading... \n\n";
            }
        }


        std::cout << "Dividing BGEN file into " << threads << " block(s)...\n";
        std::cout << "Identifying start position of each block...\n";
        std::vector<uint> endIndex(threads);
        int nBlocks = ceil(nSNPS / threads);
        uint index = 0;
        uint k = 0;
        uint sucessCount = 0;
        Mbgen_begin.resize(threads);
        Mbgen_end.resize(threads);
        bgenVariantPos.resize(threads);
        keepVariants.resize(threads);

        for (uint t = 0; t < threads; t++) 
        {
            endIndex[t] = ((t + 1) == threads) ? nSNPS - 1 : floor(((nSNPS / threads) * (t + 1)) - 1);
        }



        FILE* fin = bgen.fin;
        fseek(fin, offset + 4, SEEK_SET);
        for (uint snploop = 0; snploop < Mbgen; snploop++) 
        {
            long long unsigned int prev = ftell(fin);

            uint Nrow;
            if (Layout == 1) {
                ret = fread(&Nrow, 4, 1, fin); 
                if (Nrow != Nbgen) 
                {
                    std::cerr << "\nERROR: Number of samples (" << Nrow << ") with genotype probabilities does not match number of samples specified in BGEN file (" << Nbgen << ").\n\n";
                    exit(1);
                }
            }

            ushort LS; 
            ret = fread(&LS, 2, 1, fin);

            ret = fread(snpID, 1, LS, fin);
            snpID[LS] = '\0';
            if (checkSNPID) {
                if ((checkInclude) && (includeVariant.find(snpID) != includeVariant.end())) 
                {
                    sucessCount++;
                    keepVariants[k].push_back(snploop);
                    if (index == (nBlocks * k)) 
                    {
                        Mbgen_begin[k] = snploop;
                        long long int curr = ftell(fin);
                        bgenVariantPos[k] = curr - (curr - (prev));
                    }
                    if (index == endIndex[k]) 
                    {
                        Mbgen_end[k] = snploop;
                        k++;
                        if (k == threads) {break;}
                    }
                    index++;
                }
            }

            ushort LR; 
            ret = fread(&LR, 2, 1, fin);

            ret = fread(rsID, 1, LR, fin); 
            rsID[LR] = '\0';
            if (checkRSID) {
                if (checkInclude && (includeVariant.find(rsID) != includeVariant.end())) 
                {
                    sucessCount++;
                    keepVariants[k].push_back(snploop);
                    if (index == (nBlocks * k)) 
                    {
                        Mbgen_begin[k] = snploop;
                        long long unsigned int curr = ftell(fin);
                        bgenVariantPos[k] = curr - (curr - (prev));
                    }
                    if (index == endIndex[k]) 
                    {
                        Mbgen_end[k] = snploop;
                        k++;
                        if (k == threads) {break;}
                    }
                    index++;
                }
            }

            ushort LC; 
            ret = fread(&LC, 2, 1, fin);

            ret = fread(chrStr, 1, LC, fin); 
            chrStr[LC] = '\0';

            uint32_t physpos; 
            ret = fread(&physpos, 4, 1, fin);

            uint16_t LKnum;
            if (Layout == 2) {
                ret = fread(&LKnum, 2, 1, fin);
                if (LKnum != 2) {
                    std::cerr << "\nERROR: " << std::string(snpID) << " is a non-bi-allelic variant with " << LKnum << " alleles. Please filter these variants for now.\n\n";
                    exit(1);
                }
            }

            uint32_t LA; 
            ret = fread(&LA, 4, 1, fin);
            ret = fread(allele1, 1, LA, fin); 
            allele1[LA] = '\0';

            uint32_t LB; 
            ret = fread(&LB, 4, 1, fin);
            ret = fread(allele0, 1, LB, fin); 
            allele0[LB] = '\0';


            if (Layout == 2) {
                if (CompressedSNPBlocks > 0) 
                {
                    uint zLen; 
                    ret = fread(&zLen, 4, 1, fin);
                    fseek(fin, 4 + zLen - 4, SEEK_CUR);

                }
                else 
                {
                    uint zLen; 
                    ret = fread(&zLen, 4, 1, fin);
                    fseek(fin, zLen, SEEK_CUR);
                }
            }
            else 
            {
                if (CompressedSNPBlocks == 1) 
                {
                    uint zLen; 
                    ret = fread(&zLen, 4, 1, fin);
                    fseek(fin, zLen, SEEK_CUR);

                }
                else 
                {
                    fseek(fin, 6 * Nbgen, SEEK_CUR);
                }
            }
        }

        if (sucessCount != nSNPS) 
        {
            std::cerr << "\nERROR: There are one or more SNPs in BGEN file with " << IDline << " not in " << includeVariantFile << ".\n\n";
            exit(1);
        }
    }
    else 
    {

        filterVariants = false;
        std::cout << "Detected " << std::thread::hardware_concurrency() << " available thread(s)...\n";
        if (Mbgen < threads) 
        {
            threads = Mbgen;
            std::cout << "Number of variants (" << Mbgen << ") is less than the number of specified threads (" << threads << ")...\n";
            std::cout << "Using " << threads << " for multithreading... \n\n";
        }
        else 
        {
            std::cout << "Using " << threads << " for multithreading... \n\n";
        }

        std::cout << "Dividing BGEN file into " << threads << " block(s)..." << std::endl;
        Mbgen_begin.resize(threads);
        Mbgen_end.resize(threads);
        bgenVariantPos.resize(threads);
        std::cout << std::flush;
        keepVariants.resize(threads);
        std::cout << std::flush;
        for (uint t = 0; t < threads-1; t++) 
        {
            Mbgen_begin[t] = floor((Mbgen / threads) * t);
            Mbgen_end[t] = floor(((Mbgen / threads) * (t + 1)) - 1);
        }
        Mbgen_begin[threads-1] = floor((Mbgen / threads) * (threads - 1));
        Mbgen_end[threads-1] = Mbgen - 1;

        uint t = 0;
        FILE* fin = bgen.fin;
        fseek(fin, offset + 4, SEEK_SET);
        for (uint snploop = 0; snploop < Mbgen; snploop++) 
        {
            if (snploop == Mbgen_begin[t]) {
                bgenVariantPos[t] = ftell(fin);
                t++;
                if (t == (Mbgen_begin.size())) 
                {
					break;
				}
            }

            uint Nrow;
            if (Layout == 1) 
            {
                ret = fread(&Nrow, 4, 1, fin);
                if (Nrow != Nbgen) 
                {
                    std::cerr << "\nERROR: Number of samples (" << Nrow << ") with genotype probabilities does not match number of samples specified in BGEN file (" << Nbgen << ").\n\n";
                    exit(1);
                }
            }

            ushort LS; 
            ret = fread(&LS, 2, 1, fin);
            ret = fread(snpID, 1, LS, fin); 
            snpID[LS] = '\0';

            ushort LR; 
            ret = fread(&LR, 2, 1, fin);
            ret = fread(rsID, 1, LR, fin); 
            rsID[LR] = '\0';

            ushort LC; 
            ret = fread(&LC, 2, 1, fin);
            ret = fread(chrStr, 1, LC, fin); 
            chrStr[LC] = '\0';

            uint32_t physpos; 
            ret = fread(&physpos, 4, 1, fin);

            uint16_t LKnum;
            if (Layout == 2) 
            {
                ret = fread(&LKnum, 2, 1, fin);
                if (LKnum != 2) 
                {
                    std::cerr << "\nERROR: " << std::string(snpID) << " is a non-bi-allelic variant with " << LKnum << " alleles. Please filter these variants for now.\n\n";
                    exit(1);
                }
            }

            uint32_t LA; 
            ret = fread(&LA, 4, 1, fin);
            ret = fread(allele1, 1, LA, fin); 
            allele1[LA] = '\0';

            uint32_t LB; 
            ret = fread(&LB, 4, 1, fin);
            ret = fread(allele0, 1, LB, fin); 
            allele0[LB] = '\0';
            // Seeks past the uncompressed genotype.
            if (Layout == 2) 
            {
                if (CompressedSNPBlocks > 0) 
                {
                    uint zLen;  
                    ret = fread(&zLen, 4, 1, fin);
                    ret = fseek(fin, 4 + zLen - 4, SEEK_CUR);

                }
                else 
                {
                    uint zLen; 
                    ret = fread(&zLen, 4, 1, fin);
                    ret = fseek(fin, zLen, SEEK_CUR);
                }
            }
            else {
                if (CompressedSNPBlocks == 1) 
                {
                    uint zLen;  
                    ret = fread(&zLen, 4, 1, fin);
                    ret = fseek(fin, zLen, SEEK_CUR);

                }
                else 
                {
                    ret = fseek(fin, 6 * Nbgen, SEEK_CUR);
                }
            }
        }        
    }
    
    std::cout << std::flush;
    (void)ret;
    delete[] snpID;
    delete[] rsID;
    delete[] chrStr;
    delete[] allele1;
    delete[] allele0;
}




/***********************************************************************************
Bgen13GetTwoVals function return probs
************************************************************************************/


void Bgen13GetTwoVals(const unsigned char* prob_start, uint32_t bit_precision, uintptr_t offset, uintptr_t* first_val_ptr, uintptr_t* second_val_ptr) {

    switch (bit_precision) {
    case 8:
        *first_val_ptr = prob_start[0];
        prob_start += offset;
        *second_val_ptr = prob_start[0];
        break;
    case 16:
        *first_val_ptr = prob_start[0] | (prob_start[1] << 8);
        prob_start += offset;
        *second_val_ptr = prob_start[0] | (prob_start[1] << 8);
        break;
    case 24:
        *first_val_ptr = prob_start[0] | (prob_start[1] << 8) | (prob_start[2] << 16);
        prob_start += offset;
        *second_val_ptr = prob_start[0] | (prob_start[1] << 8) | (prob_start[2] << 16);
        break;
    case 32:
        *first_val_ptr = prob_start[0] | (prob_start[1] << 8) | (prob_start[2] << 16) | (prob_start[3] << 24);
        prob_start += offset;
        *second_val_ptr = prob_start[0] | (prob_start[1] << 8) | (prob_start[2] << 16) | (prob_start[3] << 24);
        break;
    }

}

/***********************************************************************************
calcDosage function return dosage of genotype
************************************************************************************/

std::vector<std::vector<float>>  calcDosage(std::string bgenFile, int stream_snps, int thread_num, double sigma2, Bgen bgen) 
{
    auto start_time = std::chrono::high_resolution_clock::now();
    std::vector<std::vector<float>> dosages;
    uint maxLA = 65536;
    char* snpID   = new char[maxLA + 1];
    char* rsID    = new char[maxLA + 1];
    char* chrStr  = new char[maxLA + 1];
    char* allele1 = new char[maxLA + 1];
    char* allele0 = new char[maxLA + 1];
    uint Nbgen  = bgen.Nbgen;
    uint Layout = bgen.Layout;
    uint CompressedSNPBlocks = bgen.CompressedSNPBlocks;

    std::string physpos_tmp;
    std::vector <uchar> zBuf;
    std::vector <uchar> shortBuf;

    std:: vector <uchar> zBuf1;
    std::vector <uint16_t> shortBuf1;
    uLongf destLen1 = 6 * Nbgen;
    if (Layout == 1) {
        if (CompressedSNPBlocks == 0) {
            zBuf1.resize(destLen1);
        }
        else {
            shortBuf1.resize(destLen1);
        }
    }

    bool filterVariants = bgen.filterVariants;
    int samSize = bgen.new_samSize;
    std::vector<long int> include_idx = bgen.include_idx;
    std::vector<uint> keepVariants = bgen.keepVariants[thread_num];
    uint snploop = bgen.Mbgen_begin[thread_num], end = bgen.Mbgen_end[thread_num];

    struct libdeflate_decompressor* decompressor = libdeflate_alloc_decompressor();

    FILE* fin3;
    fin3 = fopen(bgenFile.c_str(), "rb");
    long long unsigned int byte = bgen.bgenVariantPos[thread_num];
    fseek(fin3, byte, SEEK_SET);

    int variant_index = 0;
    int keepIndex = 0;
    int ret;
    while (snploop <= end) 
    {
        int stream_i = 0;
        std::vector<float> dosageList;  
        
        while (stream_i < stream_snps) //This version only suport stream_snps==1
        {
            if (snploop == (end + 1) && stream_i == 0) 
            {
                break;
            }
            if (snploop == end + 1 && stream_i != 0) 
            {
                stream_snps = stream_i;
                break;
            }
            snploop++;


            uint Nrow;
            if (Layout == 1) 
            {
                ret = fread(&Nrow, 4, 1, fin3);
                if (Nrow != Nbgen) {
                    std::cerr << "\nERROR: Number of samples (" << Nrow << ") with genotype probabilities does not match number of samples specified in BGEN file (" << Nbgen << ").\n\n";
                    exit(1);
                }
            }

            ushort LS; 
            ret = fread(&LS, 2, 1, fin3);
            ret = fread(snpID, 1, LS, fin3); 
            snpID[LS] = '\0';

            ushort LR; 
            ret = fread(&LR, 2, 1, fin3);
            ret = fread(rsID, 1, LR, fin3); 
            rsID[LR] = '\0';

            ushort LC; 
            ret = fread(&LC, 2, 1, fin3);
            ret = fread(chrStr, 1, LC, fin3); 
            chrStr[LC] = '\0';

            uint32_t physpos; 
            ret = fread(&physpos, 4, 1, fin3);
            physpos_tmp = std::to_string(physpos);

            uint16_t LKnum;
            if (Layout == 2) 
            {
                ret = fread(&LKnum, 2, 1, fin3);
                if (LKnum != 2) {
                    std::cout << "\nERROR: " << snpID << " is a non-bi-allelic variant with " << LKnum << " alleles. Please filter these variant for now. \n\n";
                    exit(1);
                }
            }

            uint32_t LA;  
            ret = fread(&LA, 4, 1, fin3);
            ret = fread(allele1, 1, LA, fin3); 
            allele1[LA] = '\0';

            uint32_t LB; 
            ret = fread(&LB, 4, 1, fin3);
            ret = fread(allele0, 1, LB, fin3); 
            allele0[LB] = '\0';

            uint nMissing = 0;
            std::vector<uint> missingIndex;

            if (Layout == 1)
             {
                uint16_t* probs_start;
                if (CompressedSNPBlocks == 1) {
                    uint zLen; 
                    ret = fread(&zLen, 4, 1, fin3);
                    zBuf1.resize(zLen);
                    ret = fread(&zBuf1[0], 1, zLen, fin3);
                    if (libdeflate_zlib_decompress(decompressor, &zBuf1[0], zLen, &shortBuf1[0], destLen1, NULL) != LIBDEFLATE_SUCCESS) {
                        std::cerr << "\nERROR: Decompressing " << snpID << " block failed with libdeflate.\n\n";
                        exit(1);
                    }
                    probs_start = &shortBuf1[0];
                }
                else {
                    ret = fread(&zBuf1[0], 1, destLen1, fin3);
                    probs_start = reinterpret_cast<uint16_t*>(&zBuf1[0]);
                }

                // read genotype probabilities
                const double scale = 1.0 / 32768;
                int idx_k = 0;   
                       
                for (uint i = 0; i < Nbgen; i++) 
                {
                    if (include_idx[idx_k] == i) {
                        double p11 = probs_start[3 * i] * scale;
                        double p10 = probs_start[3 * i + 1] * scale;
                        double p00 = probs_start[3 * i + 2] * scale;

                        if (p11 == 0 && p10 == 0 && p00 == 0) {
                            missingIndex.push_back(idx_k);
                            nMissing++;
                        }
                        else {
                            double pTot = p11 + p10 + p00;
                            double dosage = (2 * p00 + p10) / pTot;
                            dosageList.push_back(dosage);
                        }

                        idx_k++;
                    }
                }

            } // end of reading genotype data when Layout = 1

            if (Layout == 2) 
            {
                uint zLen; 
                ret = fread(&zLen, 4, 1, fin3);

                if (filterVariants && keepVariants[keepIndex] + 1 != snploop) {
                    CompressedSNPBlocks > 0 ? fseek(fin3, 4 + zLen - 4, SEEK_CUR) : fseek(fin3, zLen, SEEK_CUR);
                    continue;
                }

                uint DLen;
                uchar* bufAt;
                if (CompressedSNPBlocks == 1) 
                {
                    zBuf.resize(zLen - 4);
                    ret = fread(&DLen, 4, 1, fin3);
                    ret = fread(&zBuf[0], 1, zLen - 4, fin3);
                    shortBuf.resize(DLen);
                    uLongf destLen = DLen;

                    if (libdeflate_zlib_decompress(decompressor, &zBuf[0], zLen - 4, &shortBuf[0], destLen, NULL) != LIBDEFLATE_SUCCESS) {
                        std::cerr << "\nERROR: Decompressing " << snpID << " block failed\n\n";
                        exit(1);
                    }
                    bufAt = &shortBuf[0];
                }
                else if (CompressedSNPBlocks == 2) 
                {
                    zBuf.resize(zLen - 4);
                    ret = fread(&DLen, 4, 1, fin3);
                    ret = fread(&zBuf[0], 1, zLen - 4, fin3);
                    shortBuf.resize(DLen);
                    uLongf destLen = DLen;

                    size_t ret = ZSTD_decompress(&shortBuf[0], destLen, &zBuf[0], zLen - 4);
                    if (ZSTD_isError(ret)) 
                    {
                        std::cout << "ZSTD ERROR: " << ZSTD_getErrorName(ret);
                    }
                    bufAt = &shortBuf[0];
                }
                else 
                {
                    zBuf.resize(zLen);
                    ret = fread(&zBuf[0], 1, zLen, fin3);
                    bufAt = &zBuf[0];
                }


                uint32_t N; 
                memcpy(&N, bufAt, sizeof(int32_t));
                if (N != Nbgen) 
                {
                    std::cerr << "\nERROR: " << snpID << " number of samples (" << N << ") with genotype probabilties does not match number of samples specified in BGEN file (" << Nbgen << ").\n\n";
                    exit(1);
                }

                uint16_t K; 
                memcpy(&K, &(bufAt[4]), sizeof(int16_t));
                if (K != 2U) 
                {
                    std::cout << "\nERROR: There are SNP(s) with more than 2 alleles (non-bi-allelic). Currently unsupported. \n\n";
                    exit(1);
                }
            
                const uint32_t min_ploidy = bufAt[6];
                if (min_ploidy != 2U) 
                {
                    std::cerr << "\nERROR: " << snpID << " has minimum ploidy " << min_ploidy << ". Currently unsupported. \n\n";
                    exit(1);
                }

                const uint32_t max_ploidy = bufAt[7];
                if (max_ploidy != 2U) 
                {
                    std::cerr << "\nERROR: " << snpID << " has maximum ploidy " << max_ploidy << ". Currently unsupported. \n\n";
                    exit(1);
                }

                const unsigned char* missing_and_ploidy_info = &(bufAt[8]);
                const unsigned char* probs_start = &(bufAt[10 + N]);
                
                const uint32_t is_phased = probs_start[-2];
                if (is_phased != 1 && is_phased != 0) 
                {
                    std::cerr << "\nERROR: " << snpID << " has phased value of " << is_phased << ". This must be 0 or 1. \n\n";
                    exit(1);
                }
                
                const uint32_t bit_precision = probs_start[-1];
                if (bit_precision != 8 && bit_precision != 16 && bit_precision != 24 && bit_precision != 32) 
                {
                    std::cerr << "\nERROR: Bits to store probabilities must be 8, 16, 24, or 32. \n\n";
                    exit(1);
                }
                const uintptr_t numer_mask = (1U << bit_precision) - 1;
                const uintptr_t probs_offset = bit_precision / 8;


                int idx_k = 0;
                if (!is_phased) 
                {
                    for (uint32_t i = 0; i < N; i++) 
                    {
                        const uint32_t missing_and_ploidy = missing_and_ploidy_info[i];
                        uintptr_t numer_aa;
                        uintptr_t numer_ab;

                        if (missing_and_ploidy == 2) 
                        {
                            Bgen13GetTwoVals(probs_start, bit_precision, probs_offset, &numer_aa, &numer_ab);
                            probs_start += (probs_offset * 2);
                        }
                        else if (missing_and_ploidy == 130) 
                        {
                            if (include_idx[idx_k] == i) 
                            {
                                nMissing++;
                                missingIndex.push_back(idx_k);
                                idx_k++;
                            }
                            probs_start += (probs_offset * 2);
                            continue;
                        }
                        else 
                        {
                            std::cerr << "\nERROR: Ploidy value " << missing_and_ploidy << " is unsupported. Must be 2 or 130 (missing). \n\n";
                            exit(1);
                        }

                        if (include_idx[idx_k] == i) {
                            double p11 = numer_aa / double(1.0 * (numer_mask));
                            double p10 = numer_ab / double(1.0 * (numer_mask));
                            double dosage = 2 * (1 - p11 - p10) + p10;
                            dosageList.push_back(dosage);
                            idx_k++;
                        }
                    }

                } 
                else 
                {
                    for (uint32_t i = 0; i < N; i++) 
                    {
                        const uint32_t missing_and_ploidy = missing_and_ploidy_info[i];
                        uintptr_t numer_aa;
                        uintptr_t numer_ab;

                        if (missing_and_ploidy == 2) 
                        {
                            Bgen13GetTwoVals(probs_start, bit_precision, probs_offset, &numer_aa, &numer_ab);
                            probs_start += (probs_offset * 2);
                        }
                        else if (missing_and_ploidy == 130) 
                        {
                            if (include_idx[idx_k] == i) {
                                nMissing++;
                                missingIndex.push_back(idx_k);
                                idx_k++;
                            }
                            probs_start += (probs_offset * 2);
                            continue;
                        }
                        else 
                        {
                            std::cerr << "\nERROR: Ploidy value " << missing_and_ploidy << " is unsupported. Must be 2 or 130. \n\n";
                            exit(1);
                        }

                        if (include_idx[idx_k] == i) 
                        {
                            double p11 = numer_aa / double(1.0 * (numer_mask));
                            double p10 = numer_ab / double(1.0 * (numer_mask));
                            double dosage = 2 - (p11 + p10);
                            dosageList.push_back(dosage);
                            idx_k++;
                        }
                    }
                }
            } // end of layout 2
          
            variant_index++;
            stream_i++;
            keepIndex++;
        } // end of stream_i

        if ((snploop == (end + 1)) & (stream_i == 0)) 
        {
            break;
        }  
        dosages.emplace_back(dosageList);
    } // end of snploop 


    libdeflate_free_decompressor(decompressor);
    delete[] snpID;
    delete[] rsID;
    delete[] chrStr;
    delete[] allele1;
    delete[] allele0;
    (void)ret;

    // Close files
    fclose(fin3);

    auto end_time = std::chrono::high_resolution_clock::now();
    std::cout << "Thread " << thread_num << " finished in ";
    // printExecutionTime1(start_time, end_time);
}

