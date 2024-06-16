import os
import subprocess
from tqdm import tqdm

# Funzione per modificare il file di input con la nuova temperatura e il nuovo campo magnetico
def modify_input_file(input_file, new_temp, h=0.0, rest=1):
    # Leggi le linee dal file di input
    with open(input_file, 'r') as f:
        lines = f.readlines()

    # Modifica le linee necessarie
    for i, line in enumerate(lines):
        if line.startswith('TEMP'):
            lines[i] = f'TEMP                   {new_temp}\n'  # Modifica la riga della temperatura
        elif line.startswith('SIMULATION_TYPE'):
            values = line.split()  # Divide la riga in una lista di stringhe
            values[-1] = str(h)    # Modifica l'ultimo elemento (campo magnetico)
            lines[i] = ' '.join(values) + '\n'  # Ricompone la riga e aggiungi un newline
        if line.startswith('RESTART'):
            lines[i] = f'RESTART                {rest}\n'  # Modifica la riga del restart

    # Scrivi le linee modificate nel file di input
    with open(input_file, 'w') as f:
        f.writelines(lines)

# Funzione per estrarre l'ultima riga dal file di output e salvarla nel file di destinazione
def extract_last_line(output_file, destination_file, new_temp):
    with open(output_file, 'r') as f:
        lines = f.readlines()
        if lines:
            last_line = lines[-1].split()  # Splitta l'ultima riga in una lista di stringhe
            last_line[0] = f'{new_temp}'   # Modifica la prima colonna (temperatura)
            last_line.pop(1)               # Rimuovi la seconda colonna
            last_line_str = ' '.join(last_line)  # Unisci la lista in una stringa
    
    # Controlla se il file di destinazione esiste
    if os.path.exists(destination_file):
        # Se il file esiste, aprilo in modalità append ('a')
        with open(destination_file, 'a') as dest:
            dest.write(last_line_str + '\n')
    else:
        # Se il file non esiste, aprilo in modalità scrittura ('w')
        with open(destination_file, 'w') as dest:
            dest.write(last_line_str + '\n')

# Funzione per cancellare il contenuto di un file
def clear_file_content(file_path):
    with open(file_path, 'w') as f:
        f.truncate(0)  # Tronca il file a lunghezza zero, eliminando il contenuto

def main():
    # Definisci le variabili
    input_file = '../INPUT/input.dat'
    cxx_program = './simulator.exe'

    # Controlla il valore di SIMULATION_TYPE all'inizio
    with open(input_file, 'r') as f:
        first_line = f.readline().strip()
        if not first_line.startswith('SIMULATION_TYPE 2'):
            print("Il parametro SIMULATION_TYPE non è impostato su 2. Modificare il file di input.")
            return

    initial_temp = 3
    num_iterations = 25

    output_dir = '../OUTPUT'
    data_dir = '../../lezione_06/OUTPUT_MNT'
    prop = ["total_energy", "specific_heat", "susceptibility"]
    
    # Pulisci il contenuto dei file di destinazione prima di scrivere
    for name in prop:
        destination_file = os.path.join(data_dir, f'{name}.dat')
        clear_file_content(destination_file)
    destination_file = os.path.join(data_dir, f'magnetization.dat')
    clear_file_content(destination_file)

    # Esegui la simulazione per le proprietà con h=0.0
    print("Running simulation for properties with h=0.0")
    with tqdm(total=num_iterations) as pbar:
        for i in range(num_iterations+1):
            new_temp = "{:.1f}".format(initial_temp - i * 0.1)
            if i == 0:
                modify_input_file(input_file, new_temp, h=0.0, rest=0)
            else:
                modify_input_file(input_file, new_temp, h=0.0, rest=1)
            subprocess.run(cxx_program)
            for name in prop:
                output_file = os.path.join(output_dir, f'{name}.dat')
                destination_file = os.path.join(data_dir, f'{name}.dat')
                extract_last_line(output_file, destination_file, new_temp)
            pbar.update(1)

    # Esegui la simulazione per la magnetizzazione con h=0.02
    h_val = 0.02
    print(f"Running simulation for magnetization with h={h_val}")
    with tqdm(total=num_iterations) as pbar:
        for i in range(num_iterations+1):
            new_temp = "{:.1f}".format(initial_temp - i * 0.1)
            if i == 0:
                modify_input_file(input_file, new_temp, h=h_val, rest=0)
            else:
                modify_input_file(input_file, new_temp, h=h_val, rest=1)
            subprocess.run(cxx_program)
            output_file = os.path.join(output_dir, f'magnetization.dat')
            destination_file = os.path.join(data_dir, f'magnetization.dat')
            extract_last_line(output_file, destination_file, new_temp)
            pbar.update(1)
    
    modify_input_file(input_file, 1, 0.0, 0)

if __name__ == "__main__":
    main()
