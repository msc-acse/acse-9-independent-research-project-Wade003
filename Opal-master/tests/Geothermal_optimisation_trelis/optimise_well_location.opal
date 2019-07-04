<?xml version='1.0' encoding='utf-8'?>
<OPAL>
  <opal_operation>
    <ga_optimisation>
      <Output_filename>
        <string_value lines="1">Optimisation_wells</string_value>
      </Output_filename>
      <Executable>
        <string_value lines="1">/home/psalinas/Documents/workspace/MultiFluids_Dev/bin/icferst</string_value>
      </Executable>
      <Input_file>
        <string_value lines="1">optimodel.mpml</string_value>
      </Input_file>
      <convergence_settings>
        <Maximum_generations>
          <integer_value rank="0">1</integer_value>
        </Maximum_generations>
        <Gradient_convergence>
          <real_value rank="0">1e-2</real_value>
        </Gradient_convergence>
      </convergence_settings>
      <Trelis_integration>
        <trelis_path>/opt/Trelis-16.3/bin/trelis</trelis_path>
        <trelis_input_file>buildmodel.jou</trelis_input_file>
      </Trelis_integration>
      <Number_processors>
        <integer_value rank="0">1</integer_value>
      </Number_processors>
      <Locations_to_study>
        <value>
          <integer_value rank="0">2</integer_value>
        </value>
        <Mind_the_gap>
          <integer_value rank="0">100</integer_value>
        </Mind_the_gap>
        <Initial_guess>100,100,900,900</Initial_guess>
      </Locations_to_study>
      <Population_generation>
        <integer_value rank="0">1</integer_value>
      </Population_generation>
      <Breeding_probability>
        <value>
          <real_value rank="0">0.5</real_value>
        </value>
      </Breeding_probability>
      <Mutation_probability>
        <value>
          <real_value rank="0">0.25</real_value>
        </value>
      </Mutation_probability>
      <Precision>
        <real_value rank="0">0.05</real_value>
      </Precision>
      <Fitness>
        <Optimisation>
          <Maximise/>
        </Optimisation>
        <producer_ids>8</producer_ids>
        <injector_ids>7</injector_ids>
        <python_function>
          <string_value lines="20" type="code" language="python">val = np.sum(prod_temp) * 4185.5 * 1000 - abs(np.sum(inject_temp) * 4185.5 * 1000)</string_value>
          <comment>maximise: water_production * temperature * Cp * density</comment>
        </python_function>
      </Fitness>
      <GA_method>
        <Evolutionary_algorithm>
          <eaMuCommaLambda>
            <Mu>
              <integer_value rank="0">8</integer_value>
            </Mu>
            <Lambda>
              <integer_value rank="0">16</integer_value>
            </Lambda>
          </eaMuCommaLambda>
        </Evolutionary_algorithm>
        <Selection_method>
          <selSPEA2/>
        </Selection_method>
        <Use_CMA>
          <centroid>
            <integer_value rank="0">20</integer_value>
          </centroid>
          <sigma>
            <integer_value rank="0">40</integer_value>
          </sigma>
        </Use_CMA>
      </GA_method>
      <Variables Variable_pattern="XXX" name="X">
        <Min_limit>
          <integer_value rank="0">100</integer_value>
        </Min_limit>
        <Max_limit>
          <integer_value rank="0">900</integer_value>
        </Max_limit>
      </Variables>
      <Variables Variable_pattern="YYY" name="Y">
        <Min_limit>
          <integer_value rank="0">100</integer_value>
        </Min_limit>
        <Max_limit>
          <integer_value rank="0">900</integer_value>
        </Max_limit>
      </Variables>
    </ga_optimisation>
  </opal_operation>
</OPAL>
