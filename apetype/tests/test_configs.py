import unittest

class test_ConfigBase(unittest.TestCase):
    def setUp(self):
        import apetype as at
        class DepDepSettings(at.ConfigBase):
            c: int
            c_b: float = 3
            
        class DepSettings(at.ConfigBase):
            depdep: DepDepSettings
            b: int
            b_b: float
            b_c: int = 4
            
        class Settings(at.ConfigBase):
            depconfig: DepSettings
            a: int
            a_b: float
            a_c: int = 4

        self.Settings = Settings

    def tearDown(self):
        del self.Settings

    def test_config_with_deps_wo_prefix(self):
        settings = self.Settings(parse=False, prefix=False)
        settings.parse_args([0,0,0,0,0])
        self.assertEqual(settings.depconfig.b, 0)
        self.assertEqual(settings.b, 0)

    def test_config_with_deps_with_prefix(self):
        settings = self.Settings(parse=False, prefix=True)
        settings.parse_args([0,0,0,0,0])
        self.assertEqual(settings.depconfig.b, 0)
        self.assertEqual(settings.depconfig_b, 0)
