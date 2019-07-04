<?xml version='1.0' encoding='utf-8'?>
<OPAL>
  <opal_operation>
    <nirom>
      <snapshots_location>
        <path>
          <string_value lines="1">snapshots</string_value>
        </path>
        <input_file>
          <string_value lines="1">twomultilat.mpml</string_value>
        </input_file>
      </snapshots_location>
      <compression_eigh>
        <field_name name="phase1::Temperature">
          <nPOD>
            <integer_value rank="0">4</integer_value>
          </nPOD>
        </field_name>
      </compression_eigh>
      <training>
        <method>GPR</method>
      </training>
    </nirom>
  </opal_operation>
</OPAL>
