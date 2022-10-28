# Why emulators?

Contributed by: Pablo Giuliani, Ruchi Garg
Edited by: Kyle Godbey

Before we dive into the Reduced Basis Method, we first must look inward, deep inside our souls, and ask a fundamentally important question: why would we even want to learn about all this?

The Reduced Basis Method, from a practitioner point of view, can be thought of as an emulator, which is an algorithm that can be trained to “mimic” or approximate expensive calculations, usually with very little loss in accuracy, but with a gigantic gain in speed. If setting up the emulator is not too much of a hassle, then in principle we should always make use of them if the accuracy loss is indeed negligible. At the end of the day, performing the same calculation 100 times faster should be a great win for almost every research group. Let us walk you through two short examples, in which the speed up gain could not only be beneficial, but decisive in the amount of science that can be produced.


## First case: calibrating and quantifying the uncertainty of a physics model

Imagine that you and your colleagues have developed a nice physics model $f(\alpha,x)$ that depends on some list of parameters $\alpha$ and explains a phenomena we observe in nature as a function of the parameter $x$. For example, you could be modeling how likely it is for a neutron to interact with a nucleus as a function of the beam energy. The parameters $\alpha$ could be modeling, for example, how strongly the neutrons interact with other particles in the nucleus.

<center><img src="Fig1WhyEmu.png" alt="Fig. 1" style="width:320px;height:450px;"></center>

![Fig1](Fig1WhyEmu.png)

Figure1: Hypothetical situation you find yourself in after starting reading this book. Your model (in red) for a given parameter set $\alpha$ is compared with experimental data (in blue) your friends from the laboratory just measured.

You would like your model to compare well with the recent experimental data your friends have just measured at the laboratory after an expensive and long campaign. For that, you will have to explore the response of your model $f(\alpha,x)$ for different values of the parameters $\alpha$ until you find one that matches the observed data well. This can take several evaluations of your model, perhaps of the order of hundreds to several thousands, a task that could become very time consuming if each evaluation takes an appreciable time to compute. Here are two hypotetical examples on opposite sides of the computational burden spectrum:

<center><img src="Fig2WhyEmu.png" alt="Fig. 2" style="width:450px;height:100px;"></center>


If instead of just finding the set of optimal parameters $\alpha$ you are interested instead in finding a Bayesian posterior distribution of them to [quantify the uncertainty](https://bandframework.github.io/) of your model, things can quickly get out of hand. 

<center><img src="Fig3WhyEmu.png" alt="Fig. 3" style="width:450px;height:350px;"></center>

Figure 2: Same as Figure 1 with the addition of an uncertainty band in orange. In this hypotetical example, such band was built by sampling 5,000,000 values of the parameter $\alpha$ and evaluating the model $f(\alpha,x)$ to compare with the experimental data in blue.

These posteriors could sometimes take millions of evaluations to compute. If you were not to run them in parallel in a computer cluster, this will turn the total computational time needed to several years in the case of the simple model, and to several millennia in the case of the complex one. In these cases, an emulator becomes almost a required resource if you would like to finish your uncertainty quantification of your model $f(\alpha,x)$ within your lifetime.

## Second case, experimental optimization and online tune-up

Now let’s imagine that you put your experimental hat and you go to the laboratory to join the heroes that carefully measure nature in all its glory. An experimental set up consists of several elements and moving parts that all have to be synchronized almost to perfection to provide reliable data to the scientific community. For example, in the case of a facility that investigates rare nuclei, one part of the experimental set up is completely dedicated to the focusing and control of the beam of particles as they traverse their way from their production, through their target where the reaction occurs, until finally arriving at the target location. Figure 3 illustrates this situation for the Separator for Capture Reactions (SECAR).

<center><img src="Fig4WhyEmu.png" alt="Fig. 4" style="width:750px;height:450px;"></center>



Figure 3: Adapted from [a recent study](https://doi.org/10.1103/PhysRevAccelBeams.25.044601) by S. A. Miskovich et al., showing the focusing and defocusing of the beam of isotopes as it passes through the dipoles, quadrupoles, and octupoles magnet configurations at SECAR. The final step of the beam path is drawn as a star on a blue target, which will be missed by the beam if the currents passing through each of the magnents in SECAR is not set correctly.

The process of finding an optimal setting for the many parameters in the apparatus controlling the beam location can be challenging, sometimes involving the use of computationally expensive simulations and the fine tune of several components by hand. For such challenges, the use of an emulator that can mimic the alignment response and reduce the effective dimensionality in the large parameter space could help tremendously in speeding up the general tunning process. Further speeding up the necessary simulation packages could enable the on-line tunning in real time of the beam location as the experimental conditions are changing during the actual measurement, maximizing the time in which useful data is being actively collected.


We hope these two examples have served as an appetizer to motivate you for further reading and studying these techniques. If so, please join us on our next chapter where we give an intuitive description on why any of these types of technique have a chance of working. In a nutshell, most of our computations and calculations involving varying parameters are highly redundant.


