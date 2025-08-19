from Bio.Seq import Seq
from Bio import SeqIO, Entrez
from urllib.error import HTTPError, URLError
import socket
import pandas as pd


# menu
def main():
    print("Welcome to DNA Mutation Stimulator!🧬")
    print("Choose an option: ")
    print("1. Fetch a DNA sequence from GenBank")
    print("2. Enter your own DNA sequence")

    # prompt the user
    while True:
        try:
            choice = int(input("Choose 1, or 2: "))
            if choice == 1:
                apply(fetching())
                break
            elif choice == 2:
                apply(own_dna())
                break
            else:
                print("❌ Invalid input.")
        except ValueError:
            print("❌ Invalid input.")


# request a GenBank accession number and show sequence details
def fetching():
    while True:
        id = input("\nEnter GenBank accession number: ").strip()
        seq_record = get_from_genbank(id)
        if seq_record:
            print("\n✅Successfully fetched.\n")
            print(f"ID: {seq_record.id}")
            print(f"Description: {seq_record.description}")
            print(f"Length: {len(seq_record.seq)}")
            return seq_record.seq
        else:
            print("❌Could not fetch record. Please try again.\n")


# fetch the sequence from GenBank
def get_from_genbank(accessions):
    try:
        Entrez.email = input("Please enter a valid email address: ")
        handle = Entrez.efetch(
            db="nucleotide", id=accessions, rettype="gb", retmode="text")
        record = SeqIO.read(handle, "genbank")
        handle.close()
        return record

    # handle possible errors
    except HTTPError as e:
        print(f"❌HTTP Error {e.code}: {e.reason}")
    except URLError as e:
        print(f"❌URL Error: {e.reason}")
    except socket.timeout:
        print("❌Request timed out.")
    except ValueError as e:
        print(f"❌Value Error: {e}")
    except Exception as e:
        print(f"❌Unexpected error: {e}")
    return None


# enter your own sequence
def own_dna():
    while True:
        try:
            sequence = input("\nPlease enter your DNA sequence: ").upper()
            if not all(base in 'ACGTWSMKRYBDHVN' for base in sequence):
                raise ValueError
            break
        except ValueError:
            print("❌Invalid sequence. Please try again.")
    return sequence


# mutation
def mutate(mutation_type, seq, pos=int, base=None):
    if mutation_type == 1:
        seq = seq[:pos] + base + seq[pos+1:]
    elif mutation_type == 2:
        seq = seq[:pos] + base + seq[pos:]
    elif mutation_type == 3:
        seq = seq[:pos] + seq[pos+1:]
    return seq


# compare proteins to display the effect of the mutation
def mutation_effect(original, mutated):
    if original == mutated:
        return "Silent mutation(no change in protein)"
    elif '*' in mutated and '*' not in original:
        return "Nonsense mutation(stop codon introduced)"
    elif mutated[0] != original[0]:
        return "Missense or frameshift mutation"
    else:
        return "Mutation changed protein"


# apply mutation
def apply(dna):
    # prompt the user for the ID of the genetic code translation table
    while True:
        try:
            table = int(input(
                "\nEnter the ID of the genetic code translation table you want to use: ").strip())
            if table == 0 or table > 33:
                raise ValueError
            break
        except ValueError:
            print("❌Invalid id.")

    protein = translate(dna, table)

    print("""\n
What kind of mutation do you want to apply?
1. Substitution
2. Insertion
3. Deletion
""", end='')
    # select the mutation type
    while True:
        try:
            mutation_type = int(input('>> ').strip())
            if mutation_type not in [1, 2, 3]:
                raise ValueError
            break
        except ValueError:
            print("❌Invalid input.")

    # select the position for mutation
    while True:
        try:
            position = int(input("\nEnter position(0-based index): ").strip())
            if position > len(dna):
                raise ValueError
            break
        except ValueError:
            print("❌Invalid input.")

    # select the new base(s)
    if mutation_type == 3:
        new_base = None
    else:
        while True:
            try:
                new_base = input(
                    "\nEnter the base(s) to add or replace: ").upper().strip()
                if not all(base in 'ACGTWSMKRYBDHVN' for base in new_base):
                    raise ValueError
                break
            except ValueError:
                print("❌Invalid input.")
    # information to report
    mutated_dna = mutate(mutation_type, dna, position, new_base)
    mutated_protein = translate(mutated_dna, table)
    effect = mutation_effect(dna, mutated_dna)
    if mutation_type == 1:
        type = 'Substitution'
    elif mutation_type == 2:
        type = 'Insertion'
    elif mutation_type == 3:
        type = 'Deletion'

    # *****REPORT*****
    print("\n---------------------------------------------------------------Report---------------------------------------------------------------")
    data = pd.DataFrame(
        {'': [f'{type}', f'{position}', f"{new_base}", f"{effect}"]},
        index=["Mutation Type", "Position",
               "Base", "Effect of Mutation"]
    )
    df = pd.DataFrame(data)
    print(df)

    print(f"""\n
>>Original DNA: \n{dna}
\n>>Mutated DNA: \n{mutated_dna}
\n>>Original Protein: \n{protein}
\n>>Mutated Protein: \n{mutated_protein}
        """)


# translation
def translate(dna: str, table_id: int = 1):
    seq = Seq(dna)
    protein = seq.translate(table=table_id)
    return str(protein)


if __name__ == "__main__":
    main()
