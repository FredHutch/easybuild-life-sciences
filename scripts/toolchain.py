
class Toolchain():
    """
        methods for working with EasyBuild toolchains
    """
    def __init__(self, tc):
        self.tc = tc
        self.tc_min = 0
        self.tc_max = 11
        self.toolchains = [
            ['fosscuda-2019a', 'foss-2019a', 'GCCcore-8.2.0'],
            ['fosscuda-2019b', 'foss-2019b', 'GCCcore-8.3.0'],
            ['fosscuda-2020a', 'foss-2020a', 'GCCcore-9.3.0', 'GCC-9.3.0', 'gompi-2020a'],
            ['fosscuda-2020b', 'foss-2020b', 'GCCcore-10.2.0', 'GCC-10.2.0', 'gompi-2020b'],
            ['fosscuda-2021a', 'foss-2021a', 'GCCcore-10.3.0', 'GCC-10.3.0', 'gompi-2021a'],
            ['foss-2021b', 'GCCcore-11.2.0', 'GCC-11.2.0', 'gompi-2021b'],
            ['foss-2022a', 'GCCcore-11.3.0', 'GCC-11.3.0', 'gompi-2022a'],
            ['foss-2022b', 'GCCcore-12.2.0', 'GCC-12.2.0', 'gompi-2022b', 'gfbf-2022b'],
            ['foss-2023a', 'GCCcore-12.3.0', 'GCC-12.3.0', 'gompi-2023a', 'gfbf-2023a'],
            ['foss-2023b', 'GCCcore-13.2.0', 'GCC-13.2.0', 'gompi-2023b', 'gfbf-2023b'],
            ['foss-2024a', 'GCCcore-13.3.0', 'GCC-13.3.0', 'gompi-2024a', 'gfbf-2024a'],
        ]
        self.tc_versions = {
            '8.2.0': 0, '2019a': 0,
            '8.3.0': 1, '2019b': 1,
            '9.3.0': 2, '2020a': 2,
            '10.2.0': 3, '2020b': 3,
            '10.3.0': 4, '2021a': 4,
            '11.2.0': 5, '2021b': 5,
            '11.3.0': 6, '2022a': 6,
            '12.2.0': 7, '2022b': 7,
            '12.3.0': 8, '2023a': 8,
            '13.2.0': 9, '2023b': 9,
            '13.3.0': 10, '2024a': 10,
        }
        if tc in self.tc_versions:
           self.tc_index = self.tc_versions[tc]

    def cutoff(self, mod_name):
        """ return true if mod_name has toolchain GE to tc
        """
        for i in range(0,self.tc_max):
            for tc_name in self.toolchains[i]:
               if tc_name in mod_name:
                 return i >= self.tc_versions[self.tc]
        return True 

    def tc_trim(self, suffix):
        """ return <suffix> with out the toolchain text
        input: '0.58.1-foss-2023a' output: '0.58.1'
        """
        for tc_name in self.toolchains[self.tc_index]:
            if tc_name in suffix:
               return suffix.removesuffix(tc_name)[:-1]

    def tc_filter(self, mod_name):
        """ it <mod_name> is in toolchain return true
        how to deal with SYSTEM modules? They should return True
          - search all toolchains if none is found return True?
        """
        for tc_name in self.toolchains[self.tc_index]:
            if tc_name in mod_name:
                return True
        return False

    def tc_lt(self, Version):
        """ test if Version is less than <tc> ToolChain
            return True False
            example tc = 2023a  Version = 'GCC-12.3.0' == false
                    tc = 2023a  Version = 'GCC-11.2.0' == True
        """
        for tc in range(0, self.tc_index):
            for tc_name in self.toolchains[tc]:
               if tc_name in Version:
                  return True
        return False

    def tc_ge(self, Version):
        """ test if Version is less or equal to <tc> ToolChain
            return True False
            example tc = 2023a  Version = 'GCC-12.3.0' == True
        """
        for tc in range(self.tc_index, self.tc_max):
            for tc_name in self.toolchains[tc]:
               if tc_name in Version:
                  return True
        return False
