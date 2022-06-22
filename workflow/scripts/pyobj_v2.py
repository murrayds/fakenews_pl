from numpy import *
import json, os
# import meshpy.geometry as mg
from scipy.interpolate import interp1d

from numpy.linalg import norm

def expand(pts, n):
    s = pts.shape
    x = linspace(0,s[0]-1,n)
    xp = arange(s[0])
    return array([interp(x,xp,i) for i in pts.T]).T
    #return array([interp1d(x, y, kind='cubic') interp(x,xp,i) for i in pts.T]).T

class tube:
    def __init__(self,pts, ngon=3, radius=1., num='auto', name='test-tube'):
        self.get_pts(pts,num)
        self.radius = array([radius]).T
        self.ngon = ngon
        self.name = name
        self.get_mesh()
        self.max_idx = 0

    def get_pts(self,pts,num):
        if num=='auto':
            self.pts = array(pts)
        else:
            self.pts = expand(array(pts),num)

    def get_mesh(self):
        # print "calculating path normals ..."
        self.get_normals()
        # print "making cross sections ..."
        self.get_crosssections(self.ngon,self.n1,self.n2)
        # print "getting faces ..."
        #self.get_polygon(ngon)
        self.get_faces()
        #self.make_obj()

    def get_normals(self):
        # get unit vector from cross
        self.dp = dp = gradient(self.pts)[0] # to endure length is still len(pts) self.pts[1:]-self.pts[:-1]
        # find vectors v0 guaranteed to have cross prod with dp
        v0 = zeros_like(dp)
        v0[list(zip(*(enumerate(argmin(dp, axis=-1)))))] = 1
        # find two set of normals, n1 & n2, to dp
        n1 = cross(dp,v0)
        n1 /= norm(n1,axis=-1,keepdims=1)+1e-6
        n2 = cross(dp,n1)
        n2 /= norm(n2,axis=-1,keepdims=1)+1e-6
        self.n1 = n1
        self.n2 = n2

    def get_crosssections(self,n,v1,v2):
        "shape (#pts, #ngon, 3)"
        self.cross_sections = array([self.radius*(cos(2*pi*i/n)*v1+ sin(2*pi*i/n)*v2) +self.pts \
                                     for i in range(n)]).transpose((1,0,2))

    def get_faces(self):
        g = self.ngon
        self.faces = {'caps':[],'tri':[]}
        self.faces['caps'] = arange(g)[newaxis,:]+1+array([[0],[(self.pts.shape[0]-1)*g]])

        p = self.cross_sections - self.pts[:,newaxis,:]
        p0 = p[:-1,0,:]
        ix = argmin(norm(p[1:,:,:]-p0[:,newaxis,:],axis=-1), axis=-1)
        for n in range(len(p)-1):
            i = ix[n]
            for j in arange(g):
                self.faces['tri'] += [[n*g+1+j, n*g+1+((1+j)%g), (n+1)*g+1+((j+i) % g) ],
                                [(n+1)*g+1+((j+i+1) % g),(n+1)*g+1+((j+i) % g), n*g+1+((1+j)%g)]]
        self.faces['tri'] = array(self.faces['tri'])

    def make_obj_vf(self):
        sv = ""
        ii = self.max_idx + 0
        for p in self.cross_sections:
            for p0 in p:
                sv += "v %.4f %.4f %.4f\n" %tuple(nan_to_num(p0))
                ii +=1
        self.sv = sv
        # end caps
        sf = ""
#         for j in self.faces['caps']+self.max_idx:
#             sf += "f " +' '.join(['%d'%(i) for i in j])+'\n'
        j1,j2 = self.faces['caps']+self.max_idx
        sf += "f " +' '.join(['%d'%(i) for i in j1[::-1]])+'\n'
        sf += "f " +' '.join(['%d'%(i) for i in j2])+'\n'

        for j in self.faces['tri']+self.max_idx:
                sf += "f %d %d %d\n" %tuple(j)
        self.sf = sf
        self.max_idx = ii


    def make_obj(self):
        g = self.ngon
        s = ""
        for p in self.cross_sections:
            for p0 in p:
                s+= "v %.4f %.4f %.4f\n" %tuple(nan_to_num(p0))
        # end caps
        s0 = "f " +' '.join(['%d'%(i) for i in arange(g)[::-1] + 1 ])+'\n'
        print( s0)
        s += s0
        s += "f " +' '.join(['%d'%(i) for i in arange(g) + (self.pts.shape[0]-1)*g +1 ]) +'\n'

        # find which point on the next cross-sect is closest
        p = self.cross_sections - self.pts[:,newaxis,:]
        p0 = p[:-1,0,:]
        ix = argmin(norm(p[1:,:,:]-p0[:,newaxis,:],axis=-1), axis=-1)
        for n in range(len(p)-1):
            i = ix[n]
            for j in arange(g):
                s += "f %d %d %d\n" %(n*g+1+j, n*g+1+((1+j)%g), (n+1)*g+1+((j+i) % g) )
                s += "f %d %d %d\n" %((n+1)*g+1+((j+i+1) % g),(n+1)*g+1+((j+i) % g), n*g+1+((1+j)%g))

        self.obj = s

    def save(self):
        print("Links filename: " + self.name)
        f = open(self.name+'.obj','w')
        f.write(self.obj)
        f.close()


# import meshpy.geometry as mg
class spheres_old:
    def __init__(self, pts, rs = 1., detail = 5):
        self.detail = detail
        self.pts = array(pts)
        self.rs = (rs if hasattr(rs, '__iter__') else [rs]*len(pts))
        self.get_nodes()
        self.max_idx = 0
        self.name = 'nodes'
    def get_nodes(self):
        self.nodes = []#{'vs':[],'fc':[]}
        for r,p in zip(self.rs,self.pts):
            ex = mg.make_ball(r , self.detail)
            pt = array(ex[0])+p
            fc = ex[1]
            self.nodes += [(pt,fc)]

    def make_obj_vf(self):
        sv,sf = "",""
        ii = self.max_idx
        for pt,fc in self.nodes:
            for v in pt:
                sv += "v %.4f %.4f %.4f\n" %tuple(nan_to_num(v))
                ii +=1
            for f in fc:
                #print f
                sf += "f " +' '.join(['%d'%(i) for i in array(f[0]) + 1 + self.max_idx ])+'\n'
            # after all fces added, update max_idx do that upcoming faces are indexed correctly
            self.max_idx = ii
        self.sv = sv
        self.sf = sf

    def save(self):
        self.obj = self.sv + self.sf
        print("Nodes filename: " + self.name)
        f = open(self.name + '.obj','w')
        f.write(self.obj)
        f.close()





def save_obj(fnam,tube_cross = 5 ,tu_rad = 1.,rad_n = 3,detail=5, segs=10, pth = './', out_name = None):
    name = (out_name if out_name else fnam.split('/')[-1].split('.json')[0])
    if type(pth)==dict:
        nam_node = pth['node']+name
        nam_link = pth['link']+name
    else:
        nam_node = pth+name
        nam_link = pth+name

    net = json.load(open(fnam,'r'))
    degs = get_degrees_RS(net)
    save_nodes_obj(net['nodes'],nam_node,net['info']['nodes']['radius']*degs,rad_n,detail=detail)
    save_links_obj(net['links'],nam_link,tube_cross,tu_rad,segs=segs)

def save_obj_old(fnam,tube_cross = 6 ,tu_rad = 1.,rad_n = 3):
    nam0 = fnam.split('.json')[0]

    net = json.load(open(fnam,'r'))
    save_nodes_obj(net['nodes'],nam0,net['info']['nodes']['radius'],rad_n)
    save_links_obj(net['links'],nam0,tube_cross,tu_rad)

def save_nodes_obj(nodes,name,r,rad_n,detail=5):
    sf = spheres([[0,0,0], [1,2,3]], rs = [1, 2.5],detail=detail)
    sf.max_idx = 0
    sf.pts = array(nodes['positions'])
    if hasattr(r,'__iter__'):
        sf.rs = array(r) * rad_n
    else:
        sf.rs = array(len(sf.pts)*[r * rad_n])
    sf.get_nodes()
    print ("making node obj...")
    sf.make_obj_vf()
    sf.name = name + '-nodes'
    sf.save()

def save_links_obj(links, name, tube_cross=6,tu_rad=1, segs=10):
    s_v, s_f = '',''
    tu = tube([[0,0,0],[1,1,1]], ngon=tube_cross, num=segs)
    tu.max_idx = 0
    nl = links
    for k in nl:
        l = nl[k]
        tu.radius = l['radius'] * tu_rad
        tu.get_pts(array(l['points']), num=segs)
        tu.get_mesh()
        tu.make_obj_vf()
        s_v += tu.sv
        s_f += tu.sf

    tu.obj = s_v + s_f
    tu.name = name + '-links'

    tu.save()
#

# get node degrees
def get_degrees_RS(net):
    nl = net['links']
    ne = {}
    for k in nl:
        p = nl[k]['end_points']
        r = nl[k]['radius']
        #print(k,p,r)
        ne.setdefault(p[0],[])
        ne.setdefault(p[1],[])
        ne[p[0]]+=[r]
        ne[p[1]]+=[r]

    rs = lambda x: sqrt(sum(array(x)**2))

    degsRS = array([rs(ne[i]) for i in sorted(ne)])
    return degsRS

from pywavefront import Wavefront

def obj_vertices_faces(fnam):
    """load obj file, using pywavefront, and extract vertices and faces."""
    scene = Wavefront(fnam,create_materials=True, collect_faces=True)
    v = array(scene.vertices)
    ms = scene.mesh_list[0]
    f = array(ms.faces)
    return v, f

ball = obj_vertices_faces('./dummy_node.obj')

# import meshpy.geometry as mg
class spheres:
    def __init__(self, pts, rs = 1., detail = 5):
        self.pts = array(pts)
        self.rs = (rs if hasattr(rs, '__iter__') else [rs]*len(pts))
        self.get_nodes()
        self.max_idx = 0
        self.name = 'nodes'

    def get_nodes(self):
        self.nodes = []#{'vs':[],'fc':[]}
        for r,p in zip(self.rs,self.pts):
            if r != 0:
                self.nodes += [self.make_ball(p,r)]

    def make_ball(self, p, r):
        v,f = ball
        pt = r*v/2 +p # initial radius is 2; move origin to p
        # print(v.shape, f.shape)
        return pt, f

    def make_obj_vf(self):
        sv,sf = "",""
        ii = self.max_idx
        for pt,fc in self.nodes:
            for v in pt:
                sv += "v %.4f %.4f %.4f\n" %tuple(nan_to_num(v))
                ii +=1
            for f in fc:
                # print (f)
                sf += "f " +' '.join(['%d'%(i) for i in array(f) + 1 + self.max_idx ])+'\n'
            # after all fces added, update max_idx do that upcoming faces are indexed correctly
            self.max_idx = ii
        self.sv = sv
        self.sf = sf

    def save(self):
        self.obj = self.sv + self.sf
        f = open(self.name +'.obj','w')
        f.write(self.obj)
        f.close()
