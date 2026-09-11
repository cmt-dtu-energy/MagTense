function model = getGeneralComsolCompareBegin()

import com.comsol.model.*
import com.comsol.model.util.*

model = ModelUtil.create('Model');
model.component.create('comp1', true);
model.component('comp1').geom.create('geom1', 3);
model.component('comp1').geom('geom1').geomRep('comsol');

model.component('comp1').mesh.create('mesh1');

if (ispc)
    ModelUtil.showProgress(true);
end

end
