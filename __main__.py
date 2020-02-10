'''

'''
import pandas
import matplotlib.pyplot as plt

def ArgsInput():
    '''

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

    '''
    meta_df_ = pandas.read_csv(_meta_csv, header = None)

    print("READ METADATA\n")
    return meta_df_


def ReadRef(_ref):
    '''

    '''
    ref_df_ = pandas.DataFrame(columns = ['AA', 'FREQ'])

    for i in range(len(_ref)):
        ref_df_.loc[i] = [_ref[i], 0]

    return ref_df_


def ReadPep(_pep_csv, _pep_name):
    '''

    '''
    pep_df_raw = pandas.read_csv(_pep_csv)

    pep_df_ = pep_df_raw[pep_df_raw['Protein name']\
    .isin([_pep_name])].reset_index(drop = 'TRUE')

    print("READ CSV: " + _pep_csv)
    return pep_df_


def AlignPep(_ref, _ref_df_, _pep_df):
    '''

    '''
    for pep in _pep_df['Peptide']:
        if pep in _ref:
            ip = _ref.index(pep)
            nump = int(_pep_df[_pep_df['Peptide'].isin([pep])]['Peptide의 개수'])

            for i in range(len(pep)):
                _ref_df_.loc[ip + i, 'FREQ'] += nump

    print("ALIGNED!")
    return _ref_df_


def Plotting(_ref_df, _pep_name, _out_path, _file_name):
    '''

    '''
    for i in range(len(_ref_df.index)):
        _ref_df.loc[i, 'AA'] = '%s(%s)' % (_ref_df.loc[i, 'AA'], i)

    plt.figure(figsize = [20, 10])

    plt.bar(_ref_df['AA'], _ref_df['FREQ'])
    plt.xticks(rotation = 90)

    plt.xlabel(_pep_name)
    plt.ylabel('Frequency')

    plt.savefig(_out_path + _file_name + '.pdf')

    print("PLOTTED!: " + _out_path + _file_name + ".pdf\n")


def main():
    '''

    '''
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
        Plotting(_ref_df, _pep_name, _out_path, _file_name)


if __name__ == "__main__":
    main()
