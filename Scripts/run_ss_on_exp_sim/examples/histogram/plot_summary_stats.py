import matplotlib.pyplot as plt
import pandas as pd


"""
take in a csv produced by get_sum_stats.py and plot histograms(distributions) of the summary stats
Henry: quite old, does not have the new sum stats
"""

def main():
    df = pd.read_csv('exp_summary_stats.csv')
    length = len(df.index)
    df = df.drop(length-1) # remove mean
    df = df.drop(length-2)
    df = df.drop(length-3) # remove std
    df = df.drop(length-4)
    print(df)

    # get the summary stats names we need
    """sum_stat = list(df.columns)
    sum_stat = sum_stat[3:]
    print(sum_stat)"""

    # plot histogram
    ax = df['Aspect Ratio'].hist(by=df['sim/exp'], edgecolor='black', bins=10)
    ax[0].set_title('simulation')
    ax[1].set_title('experimental')
    plt.xlabel('Aspect Ratio')
    plt.savefig('plots/Aspect Ratio')
    plt.show()
    

    ax = df['Anisotropy'].hist(by=df['sim/exp'], edgecolor='black', bins=10)
    ax[0].set_title('simulation')
    ax[1].set_title('experimental')
    plt.xlabel('Anisotropy')
    plt.savefig('plots/Anisotropy')
    plt.show()

    ax = df['Density'].hist(by=df['sim/exp'], edgecolor='black', bins=10)
    ax[0].set_title('simulation')
    ax[1].set_title('experimental')
    plt.xlabel('Density')
    plt.savefig('plots/Density')
    plt.show()


    """df.groupby(['sim/exp'])['Aspect Ratio'].plot(kind='kde')
    plt.legend(['sim', 'exp'])
    plt.xlabel('Aspect Ratio')
    plt.show()

    df.groupby(['sim/exp'])['Aspect Ratio'].plot(kind='hist', edgecolor='black')
    plt.legend(['sim', 'exp'])
    plt.xlabel('Aspect Ratio')
    plt.show()

    df.groupby(['sim/exp'])['Anisotropy'].plot(kind='hist', edgecolor='black')
    plt.legend(['sim', 'exp'])
    plt.xlabel('Anisotropy')
    plt.show()
    
    df.groupby(['sim/exp'])['Density'].plot(kind='hist', edgecolor='black')
    plt.legend(['sim', 'exp'])
    plt.xlabel('Density')
    plt.show()"""



 

if __name__ == '__main__':
    main()
