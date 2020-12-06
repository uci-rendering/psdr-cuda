import enoki as ek

class Adam():
    """
    Implements the Adam optimizer presented in the paper *Adam: A Method for
    Stochastic Optimization* by Kingman and Ba, ICLR 2015.
    """
    def __init__(self, bsdf_map, mesh_map, bsdf_ad_keys, mesh_ad_keys,\
        lr, beta_1=0.9, beta_2=0.999, epsilon=1e-8):
        from enoki.cuda_autodiff import Float32 as Float
        # Ensure that the JIT compiler does merge 'lr' into the PTX code
        # (this would trigger a recompile every time it is changed)
        self.lr = lr 
        self.lr_v = ek.detach(Float(lr, literal=False))

        self.bsdf_map = bsdf_map
        self.mesh_map = mesh_map 
        self.bsdf_ad_keys = bsdf_ad_keys
        self.mesh_ad_keys = mesh_ad_keys
        self.beta_1 = beta_1
        self.beta_2 = beta_2
        self.epsilon = epsilon
        self.t = 0
        self.state = {}
        for k in bsdf_ad_keys:
            ek.set_requires_gradient(bsdf_map[k].reflectance.data)
            size = ek.slices(bsdf_map[k].reflectance.data)
            self.state[k] = (ek.detach(type(bsdf_map[k].reflectance.data).zero(size)),
                             ek.detach(type(bsdf_map[k].reflectance.data).zero(size)))
        for k in mesh_ad_keys:
            ek.set_requires_gradient(mesh_map[k].vertex_positions)
            size = ek.slices(mesh_map[k].vertex_positions)
            self.state[k] = (ek.detach(type(mesh_map[k].vertex_positions).zero(size)),
                             ek.detach(type(mesh_map[k].vertex_positions).zero(size)))
    def step(self):
        """ Take a gradient step """
        self.t += 1
        from enoki.cuda_autodiff import Float32 as Float
        lr_t = ek.detach(Float(self.lr * ek.sqrt(1 - self.beta_2**self.t) /
                               (1 - self.beta_1**self.t), literal=False))

        for k in self.bsdf_ad_keys:
            g_p = ek.gradient(self.bsdf_map[k].reflectance.data)
            size = ek.slices(g_p)
            assert(size == ek.slices(self.state[k][0]))
            m_tp, v_tp = self.state[k]
            m_t = self.beta_1 * m_tp + (1 - self.beta_1) * g_p
            v_t = self.beta_2 * v_tp + (1 - self.beta_2) * ek.sqr(g_p)
            self.state[k] = (m_t, v_t)
            u = ek.detach(self.bsdf_map[k].reflectance.data) - lr_t * m_t / (ek.sqrt(v_t) + self.epsilon)
            u = type(self.bsdf_map[k].reflectance.data)(u)
            ek.set_requires_gradient(u)
            self.bsdf_map[k].reflectance.data = u

        for k in self.mesh_ad_keys:
            g_p = ek.gradient(self.mesh_map[k].vertex_positions)
            size = ek.slices(g_p)
            assert(size == ek.slices(self.state[k][0]))
            m_tp, v_tp = self.state[k]
            m_t = self.beta_1 * m_tp + (1 - self.beta_1) * g_p
            v_t = self.beta_2 * v_tp + (1 - self.beta_2) * ek.sqr(g_p)
            self.state[k] = (m_t, v_t)
            u = ek.detach(self.mesh_map[k].vertex_positions) - lr_t * m_t / (ek.sqrt(v_t) + self.epsilon)
            u = type(self.mesh_map[k].vertex_positions)(u)
            ek.set_requires_gradient(u)
            self.mesh_map[k].vertex_positions = u

        
