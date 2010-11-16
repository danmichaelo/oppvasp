
from xml.dom.ext.reader import Sax2
from xml import xpath

class vasprunParser:

    def __init__(self,vasprun='vasprun.xml'):

        self.readXml(vasprun)

    def readXml(self,filename):

        doc = Sax2.FromXmlFile(filename).documentElement
        results = xpath.Evaluate( "/modeling/calculation/energy/i[@name='e_fr_energy']", doc)
        if results:
            self.toten = float(results[0].firstChild.nodeValue)
        
        print '%.6f' % (self.toten)
