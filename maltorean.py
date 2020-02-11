'''
MALDI-TOF RESULT ANALYZER by physcopatens12
===========================================

Maltorean is desinged to analyze MALDI-TOF result file from SNU Proteomics Core Facility.
Counts for each peptide read are aligned with reference protein sequences,
and converted into amino acid indexed bar graph.
Please read README.md or github.com/physcopatens12/maltorean.
'''

import pandas
import matplotlib.pyplot as plt

def ArgsInput():
    '''
    Gets arguments for metadata sheet, MALDI-TOF result sheet and output folder path.
    Input argument -h for more help.
    '''
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument('-M', '--meta', type = str, required = True,
                        help = '.csv file including metadata')
    parser.add_argument('-R', '--read', type = str, default = './',
                        help = 'folder path of read MALDI-TOF .csv files')
    parser.add_argument('-W', '--write', type = str, default = './',
                        help = 'folder path to write files')

    args_ = parser.parse_args()

    return args_


def ReadMeta(_meta_csv):
    '''
    Reads metadata sheet .csv file into pandas dataframe.
    File must be constructed with columns (without header and index):
    1) MALDI-TOF result file name (without .csv),
    2) Target reference protein name (same in result file),
    3) Amino acid sequence of reference protein.
    '''
    meta_df_ = pandas.read_csv(_meta_csv, header = None)

    print("READ METADATA\n")
    return meta_df_


def ReadRef(_ref):
    '''
    Constructs pandas dataframe of reference protein amino acid and initializes.
    '''
    ref_df_ = pandas.DataFrame(columns = ['AA', 'FREQ'])

    for i in range(len(_ref)):
        ref_df_.loc[i] = [_ref[i], 0]

    return ref_df_


def ReadPep(_pep_csv, _pep_name):
    '''
    Reads MALDI-TOF result .csv file into pandas dataframe.
    Sheets file is from Unipeptide sheet of .xlsx file.
    Columns are Protein name, Types of peptide (count) and Sum of peptide spectrum.
    '''
    pep_df_raw = pandas.read_csv(_pep_csv)

    pep_df_ = pep_df_raw[pep_df_raw['Protein name']\
    .isin([_pep_name])].reset_index(drop = 'TRUE')

    print("READ CSV: " + _pep_csv)
    return pep_df_


def AlignPep(_ref, _ref_df_, _pep_df):
    '''
    Aligns peptide sequences to reference protein.
    Each amino acid frequencies are counted by reference protein index.
    '''
    for pep in _pep_df['Peptide']:
        if pep in _ref:
            ip = _ref.index(pep)
            nump = int(_pep_df[_pep_df['Peptide'].isin([pep])]['Peptide의 개수'])

            for i in range(len(pep)):
                _ref_df_.loc[ip + i, 'FREQ'] += nump

    print("ALIGNED")
    return _ref_df_


def WriteAlign(_ref_df, _file_name, _out_path):
    '''
    Writes a sheet of reference protein amino acid and its frequency.
    '''
    _ref_df.to_csv(_out_path + _file_name + "_result.csv", index = True)

    print("WROTE: " + _out_path + _file_name + "_result.csv\n")
    return None


def Plotting(_ref_df, _pep_name, _out_path, _file_name):
    '''
    Draws bar plot of amino acid frequency indexed by reference protein sequence.
    Output size is 2000px X 5 px and formats are .pdf and .png.
    '''
    xi_aa = [0]
    for i in range(len(_ref_df.index) // 10):
        xi_aa.append((i + 1) * 10)

    plt.figure(figsize = [20, 5])

    plt.bar(_ref_df.index, _ref_df['FREQ'], color = 'black')
    plt.xticks(xi_aa, rotation = 45, fontsize = 5)

    plt.title(_pep_name + " MALDI-TOF FREQUENCY by Maltorean")
    plt.xlabel(_pep_name)
    plt.ylabel('Frequency')

    plt.savefig(_out_path + _file_name + '_result.pdf')
    plt.savefig(_out_path + _file_name + '_result.png')

    print("PLOTTED: " + _out_path + _file_name + "_result.pdf or .png")
    return None


def main():
    print ("MALDI-TOF RESULT ANALYZER by physcopatens12")
    _args = ArgsInput()
    _meta_df = ReadMeta(_args.meta)

    for i in range(len(_meta_df.index)):
        _file_name = _meta_df.iloc[i, 0]
        _pep_csv = _args.read + _file_name + '.csv'
        _pep_name = _meta_df.iloc[i, 1]
        _ref = _meta_df.iloc[i, 2]
        _out_path = _args.write

        _ref_df = ReadRef(_ref)
        _pep_df = ReadPep(_pep_csv, _pep_name)
        _ref_df = AlignPep(_ref, _ref_df, _pep_df)
        WriteAlign(_ref_df, _file_name, _out_path)
        Plotting(_ref_df, _pep_name, _out_path, _file_name)


if __name__ == "__main__":
    main()
