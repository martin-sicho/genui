import key from 'weak-key';
import React from 'react';

export function resolve(path, obj, separator='.') {
  let properties = Array.isArray(path) ? path : path.split(separator);
  return properties.reduce((prev, curr) => prev && prev[curr], obj);
}

export function TableHeaderFromItems (props) {
  const items = props.items;

  return (
    <thead>
    <tr>
      {items.map(item => (<th key={item}>{item}</th>))}
    </tr>
    </thead>
  )
}

export function TableDataFromItems (props) {
  const items = props.items;
  const dataProps = props.dataProps;
  const rowHeaderProp = props.rowHeaderProp;

  return (
    <tbody>
    {items.map(item => (
      <tr key={key(item)}>
        <th scope="row">{resolve(rowHeaderProp, item)}</th>
        {
          dataProps.map(dprop => {
            const resolved = resolve(dprop, item);
            return (
              <td key={key({dummy : resolved})}>
                {resolved}
              </td>
            )
          })
        }
      </tr>
    ))}
    </tbody>
  )
}