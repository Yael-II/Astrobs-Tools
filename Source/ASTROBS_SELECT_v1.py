import SCOPE_v2_2
import EQUATOR_v1_1
import HPMS_v1
import NGC_v1

QUIT = ["QUIT", "EXIT", "Q"]
DONE = ["DONE", "OK"]
SCOPE = ["SCOPE", "S", "1"]
EQUATOR = ["EQUATOR", "E", "2"]
HPMS = ["HPMS", "H", "3"]
NGC = ["NGC", "N", "4"]

stop = False
first = True



while not stop:
    if first:
        print("\033[35m" 
              + "{:^80}".format("=== [ ASTROBS TOOLS ] ===")
              +"\033[0m")
        first = False
    print("\033[32m"
          + "Please select a tool:"
          + "\033[0m")
    print("\033[32m"
          + "\tSCOPE: "
          + "\033[0m"
          + "Sky Coordinates for Observations Python Estimator")
    print("\033[32m"
          + "\tEQUATOR: "
          + "\033[0m"
          + "Equator Queries simbAd to create Tables of Objects")
    print("\033[32m"
          + "\tHPMS: "
          + "\033[0m"
          + "High Proper Motion Stars (Simbad Extractor)")
    print("\033[32m"
          + "\tNGC: "
          + "\033[0m"
          + "New General Catalogue (VizieR Extractor)")
    choice = input("\033[32m" + "Choice: " + "\033[0m")

    if choice.upper() in DONE + ["", " "]:
        None
    elif choice.upper() in SCOPE:
        print("\033[35m" 
              + "{:^80}".format("=== [ SCOPE ] ===")
              +"\033[0m")
        first = True
        SCOPE_v2_2.main()
    elif choice.upper() in EQUATOR: 
        print("\033[35m" 
              + "{:^80}".format("=== [ EQUATOR ] ===")
              +"\033[0m")
        first = True
        EQUATOR_v1_1.main()
    elif choice.upper() in HPMS:
        print("\033[35m" 
              + "{:^80}".format("=== [ HPMS (Simbad Extractor) ] ===")
              +"\033[0m")
        first = True
        HPMS_v1.main()
    elif choice.upper() in NGC:
        print("\033[35m" 
              + "{:^80}".format("=== [ NGC (VizieR Extractor) ] ===")
              +"\033[0m")
        first = True
        NGC_v1.main()
    elif choice.upper() in QUIT:
        stop = True
    else:
        print("\033[93m"
              + "Warning: {} is not recognized".format(choice)
              + "\033[0m")

