"use strict";

const str = "(TCAS.longest_only:0.168434,(((((BMOR.longest_only:0.155881,TRNI.longest_only:0.110983)0.203508:0.0187452,(GMEL.longest_only:0.144427,HKAH.longest_only:0.142632)0.163786:0.0266055)0.146505:0.0215354,((PPOL.longest_only:0.0747484,PXUT.longest_only:0.0455772)0.885736:0.0921985,DPLE.longest_only:0.14503)0.207377:0.0184584)0.962084:0.24131,((CNAS.longest_only:0.349197,(AGAM.longest_only:0.234182,AALB.longest_only:0.190907)0.936807:0.130355)0.376838:0.043742,(SLEB.longest_only:0.137765,(((DGRI.longest_only:0.0883319,((DHYD.longest_only:0.0388223,(DNAV.longest_only:0.0231368,(DMOJ.longest_only:0.0125944,DARI.longest_only:0.0190407)0.736652:0.0134994)0.83286:0.0288474)0.796492:0.0347524,(DVIR.longest_only:0.014752,DNOV.longest_only:0.0159103)0.921589:0.0482371)0.603559:0.0211033)0.666495:0.0256418,DBUS.longest_only:0.116091)0.626515:0.0273454,(DWIL.longest_only:0.119933,(((DOBS.longest_only:0.0349981,DGUA.longest_only:0.0393005)0.460665:0.0127046,((DPSE.longest_only:0.00592421,DPER.longest_only:0.00700677)0.568997:0.0074708,DMIR.longest_only:0.01005)0.863038:0.0308132)0.907661:0.0595581,((DBIP.longest_only:0.0312709,DANA.longest_only:0.0298577)0.914367:0.0576127,((DSER.longest_only:0.0340544,DKIK.longest_only:0.0300166)0.893216:0.0463971,((DELE.longest_only:0.037674,DRHO.longest_only:0.038239)0.576993:0.0141474,((DTAK.longest_only:0.0373668,(DSUZ.longest_only:0.0404403,DBIA.longest_only:0.0257102)0.634769:0.0168144)0.351045:0.00990442,(DEUG.longest_only:0.0422272,((DERE.longest_only:0.0220941,DYAK.longest_only:0.0249583)0.475625:0.00891292,(DMEL.longest_only:0.0160267,(DSIM.longest_only:0.0102568,DSEC.longest_only:0.00814796)0.701058:0.0113251)0.726335:0.0176944)0.755481:0.0242619)0.179262:0.00726665)0.218726:0.00853083)0.669848:0.0227083)0.461439:0.0130698)0.649987:0.0264427)0.576993:0.0245285)0.404436:0.0189147)0.288109:0.022505)0.960794:0.219423)0.796234:0.0989189)0.38406:0.0421772,((MSAC.longest_only:0.444423,BTAB.longest_only:0.375538)0.690224:0.071426,(((BTER.longest_only:0.0790374,AMEL.longest_only:0.0771587)0.931132:0.090688,(OBIR.longest_only:0.108979,(MPHA.longest_only:0.0820832,ACEP.longest_only:0.0863081)0.776631:0.0430015)0.842146:0.0705275)0.812484:0.0651028,(CSOL.longest_only:0.157669,NVIT.longest_only:0.136457)0.91643:0.099363)0.941192:0.155365)0.338406:0.0385828)1:0.168434);";

let newStr = "";

let drop = false;

for (let i = 0; i < str.length; i++) {
  if (str[i] === "." || !isNaN(parseInt(str[i]))) {
    drop = true;
  }

  if (str[i] === "," || str[i] === ")") {
    drop = false;
  }

  if (drop) {
    continue;
  }

  newStr += str[i];
}

console.log(newStr);
