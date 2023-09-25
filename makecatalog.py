from pcpy.query import makeCatalog

answer = input("\nGenerate a new catalog file (Y/N)?\n>>> ")
if answer.lower()[0] == 'y':
    makeCatalog()
    print("\nDone!")
