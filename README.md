Users to obtain read/write access to both directories:
STRI-Gamboa-ForestLandscapesFull has read/write access to both directories
garciamil@si.edu
solanom@si.edu
hernandezma3@si.edu
vasquezv@si.edu
gomezlf@si.edu

Users to obtain read access to LandscapeRaw, and read/write access to LandscapeProducts:
STRI-Gamboa-ForestLandscapesUsers has read access to the Raw directory and read/write access to the products directory.  
mullerh@si.edu
williamsoncg@si.edu
cushmank@si.edu
mitchellme@si.edu
ramosp@si.edu
cannonpg@si.edu


General conventions for ./LandscapeRaw/drone/year
The folder contains all flight missions for the Forest carbon lab. subfolder with the year are mantained for easier classification.
Each flight mission within each year follows the following convention SITE/plot/year/month/day/DRONE/mission.
		1. Years are fixed at 4 digits, months 2 digits and days 2 digits.
		2. Mission is optional for fligths covering same are for different purposes.
		3. Day corresponds to the last day of data collection. if the flights took place between september 3 and september 7. the folder will have the day 07
For example:
- BCI_50ha_2023_02_28_P4P follows the convention. 
- BCI_whole_2023_03_18_EBEE_jacaranda follows the convention

Inside each of mission folder there is 2 subdirectories that are created automatically by the server. 
	- Images subdirectory with capital I and plural -s suffix contains the images pertaining to the actual mission. 
	- The images directory can be subdivided into multiple dates, each contianing only the pictures from that date for example:
			1. 2022_03_26
			2. 2023_03_27
	- The folder Images_extra contains all the pictures non-related to the flight actual mission for example panoramics or test images.

ForestLandscapes file structure:

.
└── LandscapeRaw/
    ├── Drone/
    │   ├── BCI_50ha/
    │   │   └── 2023/
    │   │       └── 2023_01_30_BCI_50ha_P4P/
    │   │           ├── Images
    │   │           ├── Images_Extra
    │   │           ├── Telemetry
    │   │           └── Zarchive
    │   ├── BCI_AVA/
    │   ├── BCI_whole/
    │   └── DroneFlights.xlsx
    └── MLS