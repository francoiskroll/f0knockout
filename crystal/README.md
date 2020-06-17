# *crystal* fish imaging

___

Please cite the preprint if you use some of the data or code! <br />
https://doi.org/10.1101/2020.06.04.133462

___

Refers to Figure 4C, Video 1, Video 2.

Files <br />
* *eyestack_original.tif*
* *headstack_original.tif* <br />

can be found at the Zenodo version of this repository:

https://doi.org/10.5281/zenodo.3898916 <br />

## two-photon head stack

Original stack from two-photon microscope is *headstack_original.tif*.

- Opened stack in Fiji, adjusted Minimum/Maximum (Image > Adjust > Brightness/Contrast...)
- Exported in AVI (JPEG compression)
- Converted to mp4 with `ffmpeg`; command:


    ffmpeg -i headstack_jpeg.avi headstack_jpeg.mp4

- Added annotations etc in Adobe Premiere Pro.

Figure 4C: image is Z maximum projection of slice 95 to 105 (file included).

##### Two-photon parameters

___
**Note**

F0V is width/length of frame in µm. Divide by x px or y px to have µm/pixel. <br />
i.e. Frame width/height = 500 µm or 1300 px. <br />
Each px = 0.385 µm/px <br />
For scale bar: probably best is to use the whole width or height <br />
Scale bar: 130px = 50 µm
___

Software: <br />
2P_AnAb_3 <br />
z_inc 2um <br />
920nm laser <br />
gPMT 3.0V <br />
10 fr avg <br />
38 pw <br />

Acquistion parameters: <br />
x V: 3.25 <br />
y V: 3.25 <br />
x FOV: 499.99 <br />
y FOV: 499.99 <br />
x px: 1300 <br />
y px: 1300 <br />
AI: 5000 <br />
AO: 500 <br />
bin: 10 <br />
TP: 0 <br />
px_offset: -30 <br />
manShift: 0 <br />
image_offset: 1000 <br />
zPos: 150 <br />
Reps: 0 <br />
Stimuli: 0 <br />
Frame: 150 <br />
zinc: 2 <br />
LaserPower (lambda deg): 38.00 <br />
ScanType: F <br />

## two-photon eye stack

Original stack from two-photon microscope is *eyestack_original.tif*.

Same steps as for head stack to import in Adobe Premiere Pro.

##### Two-photon parameters

___
**Note**

F0V is width/length of frame in µm. Divide by x px or y px to have µm/pixel. <br />
i.e. Frame width/height = 308.24 µm or 800 px. <br />
Each px = 0.385 µm/px <br />
or each µm = 2.593 px/µm <br />
Scale bar: 129.76 px = 50 µm
___

Software: <br />
2P_AnAb_3 <br />
z_inc 2um <br />
920nm laser <br />
gPMT 3.0V <br />
10 fr avg <br />
36.5 pw <br />

Acquistion parameters: <br />
x V: 3.25 <br />
y V: 3.25 <br />
x FOV: 308.24 <br />
y FOV: 308.24 <br />
x px: 800 <br />
y px: 800 <br />
AI: 5000 <br />
AO: 500 <br />
bin: 10 <br />
TP: 0 <br />
px_offset: -30 <br />
manShift: 0 <br />
image_offset: 1000 <br />
zPos: 115 <br />
Reps: 0 <br />
Stimuli: 0 <br />
Frame: 115 <br />
zinc: 2 <br />
LaserPower (lambda deg): 36.50 <br />
ScanType: F <br />

---

Feel free to get in touch for questions

  * [![alt text][1.2]][1] [@francois_kroll](https://twitter.com/francois_kroll)

  * :email: francois@kroll.be

<!-- icons with padding -->
[1.1]: http://i.imgur.com/tXSoThF.png (twitter icon with padding)

<!-- icons without padding -->
[1.2]: http://i.imgur.com/wWzX9uB.png (twitter icon without padding)

<!-- links to your social media accounts -->
[1]: https://twitter.com/francois_kroll
