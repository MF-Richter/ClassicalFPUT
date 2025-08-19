function overprint(str)
    "This function allows to print over the last line again, i.e. no new line is created.
    -> usefull e.g. if one wants to progress through loop without unneseccarily filling the
    terminal" 
    print("\u1b[1F")
    #Moves cursor to beginning of the line n (default 1) lines up   
    print(str)   #prints the new line
   print("\u1b[0K") 
   # clears  part of the line.
   #If n is 0 (or missing), clear from cursor to the end of the line. 
   #If n is 1, clear from cursor to beginning of the line. 
   #If n is 2, clear entire line. 
   #Cursor position does not change. 

    println() #prints a new line, i really don't like this arcane codes
end