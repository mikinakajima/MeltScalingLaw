# Melt Scaling law 

This is a scaling law described in Nakajima et al. in review. To use the scaling law, please take the following steps.

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.

### Prerequisites

Make sure you have python3 on your machine and the following libraries.

```
matplotlib
numpy
sys
scipy

```

Download color maps from Fabio Crameri's website by typing

```
./DownloadColormaps.sh 
```


### Running the script

To run a script, simply type

```
python  melt_model_Nakajima_et_al.py
```
This will generate the following output

```
magma ocean depth and pressure for a melt pool model: 2278.4837243704346 km, 86.15126263895684 GPa
magma ocean depth and pressure for a global magma ocean model: 898.1255099313564 km, 29.474287387562157 GPa
magma ocean depth and pressure for a conventional model: 303.5506251133121 km, 9.21756896904488 GPa
```


```
Mtotal: Total mass normalized by Mars mass <br />
gamma: Imapctor-to-total mass ratio <br />
vel: Imapct velocity normalized by mutual escape velocity <br />
entropy0: Initial mantle entropy (before impact) <br />
angle: Impact angle (0 is a head-on collision and 90 is the most oblique impact. Choose from 0, 30, 60, 90 degrees) <br />
outputfigurename: Name of the output figure name <br />
```

And repeat

```
until finished
```

End with an example of getting some data out of the system or using it for a little demo

## Running the tests

Explain how to run the automated tests for this system

### Break down into end to end tests

Explain what these tests test and why

```
Give an example
```

### And coding style tests

Explain what these tests test and why

```
Give an example
```

## Deployment

Add additional notes about how to deploy this on a live system

## Built With

* [Dropwizard](http://www.dropwizard.io/1.0.2/docs/) - The web framework used
* [Maven](https://maven.apache.org/) - Dependency Management
* [ROME](https://rometools.github.io/rome/) - Used to generate RSS Feeds

## Contributing

Please read [CONTRIBUTING.md](https://gist.github.com/PurpleBooth/b24679402957c63ec426) for details on our code of conduct, and the process for submitting pull requests to us.

## Versioning

We use [SemVer](http://semver.org/) for versioning. For the versions available, see the [tags on this repository](https://github.com/your/project/tags). 

## Authors

* **Billie Thompson** - *Initial work* - [PurpleBooth](https://github.com/PurpleBooth)

See also the list of [contributors](https://github.com/your/project/contributors) who participated in this project.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Hat tip to anyone whose code was used
* Inspiration
* etc



# Melt Scaling Law by Nakajima et al.
Melt Scaling law from Nakajima et al. in rev.

This is a scaling law described in Nakajima et al. in review. To use the scaling law, please take the following steps.

**1. Download color maps from Fabio Crameri's website by typing**  <br />
./DownloadColormaps.sh 

2. To change input parameters, open input.txt file <br />

Mtotal: Total mass normalized by Mars mass <br />
gamma: Imapctor-to-total mass ratio <br />
vel: Imapct velocity normalized by mutual escape velocity <br />
entropy0: Initial mantle entropy (before impact) <br />
angle: Impact angle (0 is a head-on collision and 90 is the most oblique impact. Choose from 0, 30, 60, 90 degrees) <br />
outputfigurename: Name of the output figure name <br />

3. To calculate the magma ocean pressures and magma ocean sizes, simply type  <br />

python melt_model_Nakajima_et_al.py  <br />

This will calculate three different pressures;  <br />
magma ocean depth and pressure for a melt pool model: x0 km,  y0 GPa <br />
magma ocean depth and pressure for a global magma ocean model: x1 km, y1 GPa <br />
magma ocean depth and pressure for a conventional model: x2, y2 GPa <br />

These are calculated depth (in km) and pressure (in GPa) at the base of a magma ocean depth for a melt pool model, global magma ocean model, and conventional model. <br />

If you have questions, please feel free to reach out to Miki Nakajima (mnakajima@rochester.edu). If you use this scaling law, please be sure to cite our paper, Nakajima et al., in review.


