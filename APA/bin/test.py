def main():
    ff = "/mnt/data8/zhouran/proj/2020-APA_bulk/Hema_mouse/bulk_tianjin/data/bam/B_1.Log.final.out"
    with open(ff) as fh:
        for line in fh.readlines():
            line = line.strip()
            if line.startswith("Uniquely mapped reads number"):
                line = line.split("\t")
                print(line[-1])


if __name__ == "__main__":
    main()
