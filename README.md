Project title: Deleterious mutations cause evolution of lifespan extension by dietary restriction

Overview: We compared the effects of dietary restriction (DR) and DR with exposure to food odour to determine whether DR response is adaptive. We compared survival, mated and unmated reproductive output, egg size, stress resistance under heatshock, and body size under 4 different treatments:

Ad libitum food (F)
Ad libitum food with food odour (FO or F+O)
Dietary restriction for 24 hours (DR)
Dietary restriction for 24 hours with food odour (DRO or DR+O)

All data were analysed using R version 4.4.1.

File: DRO_LS
Description: Nematode survival data following treatments as described above
Variables:
  Treatment: F, F+O, DR, DR+O treatments as described in overview
  Plate.ID: 10 individuals were kept together per plate; Plate.ID represents the cohort kept       together
  Worm: Arbitrary number (1-10, of 10 worms on plate) given to individuals upon death
  Death.Date: Date that individual died
  Cause: Cause of death - W = Walled, L = Lost, D = Natural death, M = Matricide, E = Explosion
  D1: Date experiment began
  Age: Number of days individuals survived
  Event: Binary variable, where 1 = natural death, 0 = loss, walling, or any 'unnatural' deaths
  Inf: Binary variables, where 1 = infected individual, 0 = individual not infected

File: DR_Lab_repro_2
Description: Age-specific reproduction data following treatments as described above
Variables:
  ID: Individual identification number 
  Treatment: F, F+O, DR, DR+O treatments as described in overview
  D1...D8: Number of offspring produced on Days 1 through 8 

File: Male_repro_DR
Description: Age-specific reproduction data after mating and following treatments as described above
Variables:
  Block: Experimental Block (we only performed 1)
  Treatment: F, F+O, DR, DR+O treatments as described in overview
  ID: Individual identification number
  D1...D10: Number of offspring produced on Days 1 through 10
  Lost: Notes on whether individuals died or were lost, N = not lost, MAT = matricide, INF = death due to infection, L = lost

File: Mated_LS
Description: Nematodes survival data following mating and treatment as described above
Variables:
  Block: Experimental block (we only performed 1)
  Treatment: F, F+O, DR, DR+O treatments as described in overview
  Repro.ID: Plate number assigned from reproduction experiment (individuals have this only if they did not survive beyond reproduction experiment)
  LS.plate: Plate number assigned when grouped from individual reproduction plates onto lifespan plates (individuals have this only if they survived beyond reproduction experiment)
  Worm: Individual ID within Treatment
  Death.Date: Date individual died
    Cause: Cause of death - W = Walled, L = Lost, D = Natural death, M = Matricide, E = Explosion
  D1: Date experiment began
  Age: Number of days individuals survived
  Event: Binary variable, where 1 = natural death, 0 = loss, walling, or any 'unnatural' deaths
  Exp: Binary variables, where 1 = exploded individual, 0 = individual not exploded
  Inf: Binary variables, where 1 = infected individual, 0 = individual not infected
  Goop: Binary variables, where 1 = infected individual with strange goop present on plate, 0 = individual does not have goop infection
  Infect: Yes/No (Y/N) factor variable whether individual was infected or not

File: Egg_size
Description: Size of eggs produced on Days 2 and 4 of adulthood from nematodes following treatments as described in overview
Variables:
ID: Individual identification number
Treatment: F, F+O, DR, DR+O as described in overview
Day: Day of adulthood eggs were measured on (Day 2 is peak reproduction for F and F+O, Day 4 is peak reproduction for DR and DR+O)
Replicate: Number of egg measured within individual output (3 eggs were measured per individual)
Egg_size: Surface area of eggs (mm^2)

File: Binary_size
Description: Body size of individuals on Days 2 (immediately after treatment) and 4 (after 2 days of feeding following treatment) of adulthood from nematodes following treatments as described in overview
Variables:
ID: Individual identification number
Treatment: F, F+O, DR, DR+O as described in overview
Day: Day of adulthood body sizewas measured on (Day 2 or 4)
Body: Surface area of body (mm^2)

File: DR_HS
Description: Survival during and following heatshock (37C) following treatments as described in overview. Mortality checks were performed hourly during the 9 hour heatshock treatment and daily beginning the day after heatshock treatment
Variables:
Plate_ID: Plate identification number: individuals were kept in groups of ten per plate within treatments
Treatment: F, F+O, DR, DR+O as described in overview
Age: Number provided to each timepoint irrespective of hour/day (i.e., Hourly mortality check data during heatshock occurs at Age 1-9, Day 1 mortality check occurs at Age 10 (this is only used to generate a survival plot that adequately shows the survival curves during heatshock, but is not used for analysis)
Time: Time/Day that mortality was observed/recorded
Event: Binary variable where 1 = natural death, 0 = death due to unnatural causes (i.e., lost individuals, walled individuals)
Day: Day at which individuals were recorded as dead (hourly heatshock mortalities are represented as a proportion of 24 hours)

File: DR_wild2
Description: Population index data based on counts from a 5ml sample of compost from a 250 ml microcosm containing starting populations of 50 individuals following a 24-hour dietary treatment in adulthood. Samples were taken weekly over a 3-week period
Variables:
Tube: Identification number for tubes/microcosms each containing a population from a single treatment
Treatment: F, F+O, DR, DR+O as described in overview
Rep: Technical replicate produced by taking a 0.10 ml sample of liquid containing suspended nematodes from compost samples, pipetted onto an agar plate to allow for counting
Day: 7, 14, or 21, representing the day the sample was taken after microcosms were inoculated with nematode populations and placed outdoors
No_worms: Number of individual nematodes counted from a technical replicate
Block: Experiment was performed in 3 3-week intervals throughout summer of 2024
Population: Individual identification number of populations across all blocks

