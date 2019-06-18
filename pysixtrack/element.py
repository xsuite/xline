from dataclasses import dataclass

class _MetaElement(type):
    def __new__(cls, clsname, bases, dct):
        description = dct.get('_description', [])
        ann={}
        dct['__annotations__']=ann
        for name,unit,desc,default in description:
            ann[name]=object
        try:
            doc = [dct['__doc__'], '\nFields:\n']
        except KeyError:
            doc = ['\nFields:\n']
        fields = [f"{field:10} [{unit+']:':5} {desc} " for field,
                  unit, desc in description]
        doc += fields
        dct['__doc__'] = "\n".join(doc)
        newclass=super(_MetaElement, cls).__new__(cls, clsname, bases, dct)
        return dataclass(newclass)


class Base(metaclass=_MetaElement):
    def get_fields(keepextra=False):
        pass

    def to_dict(self,keepextra=False):
        return {kk: getattr(self,kk) for kk in self.get_fields(keepextra)}

    @classmethod
    def from_dict(self,keepextra=True):
        pass

    def copy(self,keepextra=True):
        return self.__class__.from_dict(self.to_dict(keepextra),keepextra)


class Element(Base):
    pass

