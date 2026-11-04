
"""
CRISPR gRNA Designer Tool
by: Ansuman Gadanayak
Scans a DNA sequence for SpCas9 NGG PAM sites, extracts the 20nt guide RNA sequence,
computes GC content, and scores/filters for optimal gRNA design.
"""

import sys
import argparse
re = __import__('re')

def calculate_gc(seq):
    """Calculate the GC percentage of a sequence."""
    seq = seq.upper()
    if not seq:
        return 0.0
    gc_count = seq.count('G') + seq.count('C')
    return round((gc_count / len(seq)) * 100, 2)

def evaluate_poly_t(seq):
    """Check for four or more consecutive thymines (Pol III terminator risk)."""
    return bool(re.search(r'TTTT', seq, re.IGNORECASE))

def scan_grnas(dna_sequence, strand="both"):
    """
    Scans sequence for NGG PAM sites (Cas9) on sense and/or antisense strands.
    """
    results = []
    dna_upper = dna_sequence.upper()
    
    def process_strand(seq, strand_label):
        for i in range(len(seq) - 22):
            if seq[i+21] == 'G' and seq[i+22] == 'G':
                full_target = seq[i:i+23] # 20nt guide + 3nt PAM
                guide = full_target[:20]
                pam = full_target[20:]
                
                if 'N' in guide:
                    continue
                
                gc = calculate_gc(guide)
                has_poly_t = evaluate_poly_t(guide)
                
                score_status = "Optimal"
                if not (40 <= gc <= 60):
                    score_status = "Suboptimal GC"
                if has_poly_t:
                    score_status = "Risk: Poly-T terminator"

                results.append({
                    "position": i + 1,
                    "strand": strand_label,
                    "guide_20nt": guide,
                    "pam": pam,
                    "full_target": full_target,
                    "gc_content": gc,
                    "status": score_status
                })

    if strand in ["plus", "both"]:
        process_strand(dna_upper, "+")
        
    if strand in ["minus", "both"]:
        trans = str.maketrans("ATCGN", "TAGCN")
        rev_comp = dna_upper.translate(trans)[::-1]
        process_strand(rev_comp, "-")
        
    return results

def main():
    print("=" * 75)
    print("           CRISPR gRNA DESIGNER TOOL - Developed by Ansuman Gadanayak")
    print("=" * 75 + "\n")

    parser = argparse.ArgumentParser(description="CRISPR SpCas9 gRNA Designer Tool")
    parser.add_argument("-i", "--input", help="Path to FASTA or text file containing the gene sequence")
    parser.add_argument("-s", "--seq", help="Raw DNA sequence string")
    parser.add_argument("--strand", choices=["plus", "minus", "both"], default="both", help="Strand to scan")
    args = parser.parse_args()

    sequence = ""
    if args.seq:
        sequence = args.seq
    elif args.input:
        try:
            with open(args.input, 'r') as f:
                lines = f.readlines()
                sequence = "".join([line.strip() for line in lines if not line.startswith(">")])
        except FileNotFoundError:
            print(f"Error: File '{args.input}' not found.")
            input("\nPress Enter to exit...")
            sys.exit(1)
    else:
        # Fallback to interactive input if no arguments are passed so the window doesn't close
        sequence = input("Enter your target DNA sequence: ").strip()

    if not sequence:
        print("Error: Empty sequence provided.")
    else:
        print(f"\n[*] Scanning sequence of length {len(sequence)} bp for SpCas9 (NGG) PAM sites...\n")
        grnas = scan_grnas(sequence, strand=args.strand)

        print(f"{'Pos':<6} {'Strand':<8} {'20nt Guide Sequence':<22} {'PAM':<6} {'GC%':<8} {'Evaluation'}")
        print("-" * 75)
        
        for g in grnas:
            print(f"{g['position']:<6} {g['strand']:<8} {g['guide_20nt']:<22} {g['pam']:<6} {g['gc_content']:<8} {g['status']}")

        print(f"\n[+] Total valid gRNA targets found: {len(grnas)}")

    # Keeps the window open until you press Enter
    input("\nPress Enter to exit...")

if __name__ == "__main__":
    main()
