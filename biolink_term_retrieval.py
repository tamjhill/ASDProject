from bmt import Toolkit
t = Toolkit()
element = t.get_element('doi')
print(element)

t.is_category('gene') 

t.is_predicate('related to') 

element_name = t.get_element_by_mapping('SEMMEDDB:CAUSES') # This returns 'causes'