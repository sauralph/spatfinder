import typer
from typing import List
from pathlib import Path
import os
import app.spa_typer as spa_typer

app = typer.Typer()

@app.command()
def main(files: List[Path], # =typer.Argument(['.'],help="files to be processed by SPA typer"),
    untypable:bool = typer.Option(True,help="include untypable patterns in the output"),
    guess:bool = typer.Option(True,help="suggest possible types for untypable patterns"),
    silent:bool = typer.Option(False,help="only produce output file, no information on stdout"),
    output:str = typer.Option('output.csv',help="name for the output file"),
    input_extension:str = typer.Option('ffn',help="file extension to include automatically when processing directory")
    ):
    """
    SPA-typer

    Detect SPA types from contigs or enssembled genomes
    """
    compile_files = []
    for path in files:
        if path.is_dir():
            compile_files.extend([ os.path.join(path, fasta) for fasta in os.listdir(path) if fasta.endswith(input_extension) ])
        if path.is_file():
            compile_files.append(path.name)
    
    # print(path.absolute)
    print(compile_files)
    spa_typer.process_files(compile_files)


if __name__ == "__main__":
    app()