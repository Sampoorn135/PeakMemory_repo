import random
import pyBigWig
from bwremote import RemoteBigWig

EX = 'https://genome.ucsc.edu/goldenPath/help/examples/bigWigExample.bw'
bw = pyBigWig.open('ex.bw')
chrom = list(bw.chroms().keys())[0]; size = bw.chroms()[chrom]
print('chrom', chrom, 'size', size)
rb = RemoteBigWig(EX)
print('remote chroms match:', rb.chroms.get(chrom)[1] == size)
print('reader initialized', flush=True)
random.seed(1)
allok = True
for _ in range(4):
    s = random.randint(0, size - 50000); e = s + random.randint(500, 5000)
    tsum = bw.stats(chrom, s, e, type='sum', exact=True)[0] or 0.0
    tcovf = bw.stats(chrom, s, e, type='coverage', exact=True)[0] or 0.0
    msum, mcov, L = rb.region_stats(chrom, s, e)
    ok = abs(tsum - msum) < 1e-2 * max(1, abs(tsum)) and abs(tcovf * L - mcov) <= 1
    allok &= ok
    print(f'{s}-{e} truth_sum={tsum:.3f} mine={msum:.3f} truth_covbp={tcovf*L:.0f} mine={mcov} OK={ok}', flush=True)
bw.close()
print('ALL OK' if allok else 'MISMATCH')
