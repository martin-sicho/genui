import React from 'react';
import { Table } from 'reactstrap';

export function PropertiesTable(props) {
  return (
    <Table size="sm" hover>
      <thead>
      <tr>
        <th>Property</th>
        <th>Value</th>
      </tr>
      </thead>
      <tbody>
      {
        props.propsList.map(propName => (
          <tr key={propName}>
            <td>{propName}</td>
            <td>{props.molWithProperties.properties[propName]}</td>
          </tr>
        ))
      }
      </tbody>
    </Table>
  )
}