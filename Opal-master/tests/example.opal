<?xml version='1.0' encoding='utf-8'?>
<Importance_mapper>
  <Output_filename>
    <string_value lines="1">Test</string_value>
  </Output_filename>
  <Executable>
    <string_value lines="1">icferst</string_value>
  </Executable>
  <Input_file>
    <string_value lines="1">test.mpml</string_value>
  </Input_file>
  <Field_to_study name="PhaseVolumeFraction">
    <field_type>
      <scalar_field/>
    </field_type>
    <perturbations_to_do>
      <integer_value rank="0">7</integer_value>
    </perturbations_to_do>
    <Sample_width>
      <real_value rank="0">0.01</real_value>
    </Sample_width>
    <Phases_to_perturb>
      <integer_value shape="2" rank="1">1 2</integer_value>
    </Phases_to_perturb>
  </Field_to_study>
  <Field_to_study name="Pressure">
    <field_type>
      <scalar_field/>
    </field_type>
    <perturbations_to_do>
      <integer_value rank="0">5</integer_value>
    </perturbations_to_do>
    <Boundary_condition/>
    <Sample_width>
      <real_value rank="0">0.01</real_value>
    </Sample_width>
    <Phases_to_perturb>
      <integer_value shape="2" rank="1">1 2</integer_value>
    </Phases_to_perturb>
  </Field_to_study>
</Importance_mapper>
