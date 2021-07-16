import sys
import requests
import json
from shutil import which
from collections import defaultdict


def post_get(fa):
    data = {
        "in_seq": fa
        }
    r = requests.post(
        'https://bio.vbio.top/api/bio/shrna-design/design/', 
        json=data
        )

    uuid = json.loads(r.text)['uuid']
    res = json.loads(
        requests.get(
            f"https://bio.vbio.top/api/bio/shrna-design/result/{uuid}/").text
        )
    return res['output']['result']


# Query_1\tENST00000659131\t100.00\t21\t0\t0\t1\t21\t1999\t1979\t4e-04\t42.1
def blast(fa, trans_id, blast, blastdb):
    import subprocess
    blastn = subprocess.Popen(
        f"{blast} -db {blastdb} -task blastn-short -outfmt 6",
        shell=True,
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT
        )
    results, err = blastn.communicate(fa.encode("utf-8"))
    results = results.decode().split('\n')

    hit_n = 0
    hit_self = 0
    if not results:
        return hit_n, hit_self
    for res in results:
        if not res:
            continue
        res = res.split('\t')
        try:
            if res[1] == trans_id:
                hit_self += 1
                continue
            if res[2] == '100.00':
                hit_n += 1
            else:
                pass
        except IndexError:
            continue
    return hit_n, hit_self


def main(fa: str, output: str, blast: str, blastdb: str):
    if which(blast) is None:
        raise ValueError(f'{blast} was not found.')
        sys.exit(0)

    fa_pool = defaultdict(str)
    cid = ''
    with open(fa) as fi, open(output, 'w') as fo:
        fo.write('\t'.join(
            ['score',
             'shRNA_seq',
             'startsite',
             'id',
             'num_of_hit',
             'num_of_hit_self']
            ) + '\n')

        for line in fi:
            if line.startswith('>'):
                cid = line.strip()[1:]
            else:
                fa_pool[cid] += line.strip()
        for label, seq_ in fa_pool.items():
            res = post_get(seq_)
            for sub_shRNA in res:
                score, seq_, st = list(sub_shRNA.values())
                hit_n, hit_self = blast(seq_, label.split('_')[0], blast, blastdb)
                fo.write('\t'.join(
                    [str(score),
                     seq_,
                     str(st),
                     label,
                     str(hit_n),
                     str(hit_self)]
                    ) + '\n')


if __name__ == '__main__':
    import fire
    fire.Fire(main)
