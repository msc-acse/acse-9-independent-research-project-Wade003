<?xml version='1.0' encoding='utf-8'?>
<OPAL>
  <opal_operation>
    <importance_map>
      <Output_filename>
        <string_value lines="1">Test</string_value>
      </Output_filename>
      <Executable>
        <string_value lines="1">~/Packages/MultiFluids_Dev/bin/multiphase_prototype</string_value>
      </Executable>
      <Input_file>
        <string_value lines="1">BL_fast.mpml</string_value>
      </Input_file>
      <functional>
        <end_time/>
        <location_of_interest>
          <real_value shape="2" rank="1">1.0 0.066666667</real_value>
        </location_of_interest>
        <type>
          <standard/>
        </type>
        <Functional_field name="Temperature">
          <field_type>
            <scalar_field/>
          </field_type>
        </Functional_field>
      </functional>
      <Field_to_study name="Velocity">
        <field_type>
          <scalar_field/>
        </field_type>
        <perturbations_to_do>
          <integer_value rank="0">6</integer_value>
        </perturbations_to_do>
        <Initial_condition/>
        <Gram-Schmidt/>
        <smoothing>
          <integer_value rank="0">6</integer_value>
        </smoothing>
        <use_G>
          <integer_value rank="0">1</integer_value>
        </use_G>
        <Sigma>
          <real_value rank="0">0.01</real_value>
        </Sigma>
        <Phases_to_perturb>
          <integer_value shape="1" rank="1">1</integer_value>
        </Phases_to_perturb>
      </Field_to_study>
    </importance_map>
  </opal_operation>
</OPAL>
